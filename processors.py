import awkward as ak
from coffea import hist, processor
import pandas as pd
import json
from dask.distributed import get_worker
from subprocess import call

# register our candidate behaviors
from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')

call("mkdir -p logs", shell=True)
#file_handler = logging.FileHandler('logs/processor.log')
file_handler = logging.FileHandler('processor.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

class VBFHHggtautauProcessor(processor.ProcessorABC):
    '''
    https://github.com/CoffeaTeam/coffea/blob/65978ea299eee653f01ccfbc11af9f0339ffc7c2/tests/test_dask_pandas.py
    https://github.com/CoffeaTeam/coffea/blob/bbfd34414530f981529ccecf16e1585d293c5389/coffea/processor/test_items/NanoTestProcessorPandas.py
    '''
    lumis = {"2016": 35.9, "2017": 41.5, "2018": 59.8}

    def __init__(self, output_path, job_tag):
        self.output_path = output_path
        self.job_tag = job_tag

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):

        output = pd.DataFrame() 
        dataset = events.metadata['dataset'] # identifier of data/MC 
        
        logger.info(f"start to process events")

        from utils import wrap_items
        item_names_to_be_wrapped = ["selectedPhoton", "tau", "electron", "muon", "jet", "tautauSVFitLoose"]
        photons, taus, electrons, muons, jets, ditau_svfit = wrap_items(events, item_names_to_be_wrapped)

        ### object selection
        from utils import select_tau, select_electron, select_muon, select_jet, select_vbf_jet

        dipho_presel = photons[:,0] + photons[:,1]
        photonID_cut = (photons.mvaID[:,0] > -0.7) & (photons.mvaID[:,1] > -0.7)
        photon_ptom_cut = (photons.pt[:,0]/dipho_presel.mass > 0.33) & (photons.pt[:,1]/dipho_presel.mass > 0.25) 
        diphoton_mass_cut = (dipho_presel.mass > 100) & (dipho_presel.mass < 180)
        extra_dipho_cut  = photonID_cut & photon_ptom_cut & diphoton_mass_cut
        events["selectedPhoton"] = photons
        events["diphoton"] = dipho_presel

        selectedMuons = select_muon(muons, photons)
        selectedElectrons = select_electron(electrons, photons)
        selectedTaus = select_tau(taus, photons, selectedElectrons, selectedMuons)
        selectedJets = select_jet(jets, photons, selectedMuons, selectedElectrons, selectedTaus) 
        selectedVBFJets = select_vbf_jet(jets, photons, selectedMuons, selectedElectrons, selectedTaus) 
        ### add selected objects back to the event ak array
        events["selectedMuon"] = selectedMuons
        events["selectedElectron"] = selectedElectrons
        events["selectedTau"] = selectedTaus
        events["selectedJet"] = selectedJets
        events["selectedVBFJet"] = selectedVBFJets

        ### event selection
        sum_charge = ak.sum(selectedElectrons.charge, axis=1) + ak.sum(selectedMuons.charge, axis=1) + ak.sum(selectedTaus.charge, axis=1)
        n_electrons = ak.num(selectedElectrons)
        n_muons = ak.num(selectedMuons)
        n_taus = ak.num(selectedTaus)
        # Require >= 1 lep/tau
        n_leptons_and_taus = n_electrons + n_muons + n_taus
        # os if 2 hadTaus or 1 hadTau + 1 lep
        os_cut = (n_leptons_and_taus == 2) & (sum_charge == 0)
        one_tau_zero_lep_cut = (n_taus == 1) & (n_electrons + n_muons == 0)
        # analysis pre-selection
        all_cuts = (extra_dipho_cut) & (os_cut | one_tau_zero_lep_cut)

        n_vbf_jets = ak.num(selectedVBFJets)

        #### save relevant info
        output["run"] = ak.to_numpy(events["run"])
        output["event"] = ak.to_numpy(events["event"])

        output["diphoton_mass"] = events["diphoton"].mass
        output["diphoton_ptom"] = events["diphoton"].pt/events["diphoton"].mass
        output["diphoton_eta"] = events["diphoton"].eta
        output["diphoton_deltaR"] = events["selectedPhoton"][:,0].delta_r(events["selectedPhoton"][:,1])
        output["lead_pho_ptom"] = events["selectedPhoton"][:,0].pt/events["diphoton"].mass
        output["lead_pho_eta"] = events["selectedPhoton"][:,0].eta
        output["lead_pho_id"] = events["selectedPhoton"][:,0].mvaID
        output["sublead_pho_ptom"] = events["selectedPhoton"][:,1].pt/events["diphoton"].mass
        output["sublead_pho_eta"] = events["selectedPhoton"][:,1].eta
        output["sublead_pho_id"] = events["selectedPhoton"][:,1].mvaID
        output["nHadTaus"] = n_taus
        output["nElectrons"] = n_electrons
        output["nMuons"] = n_muons

        from utils import getcosthetastar_cs, helicityCosTheta
        output["costhetastar_cs"] = getcosthetastar_cs(dipho_presel, ditau_svfit) 
        output["svfit_mass"] = ditau_svfit.mass
        output["costheta_gg"] = helicityCosTheta(dipho_presel, photons[:,0])
        output["costheta_tautau"] = -999
        
        cat_tautau = (n_taus == 2)
        tau1_cat_tautau = selectedTaus[cat_tautau][:,0] 
        tau2_cat_tautau = selectedTaus[cat_tautau][:,1] 
        ditau_cat_tautau = tau1_cat_tautau + tau2_cat_tautau
        output.loc[ak.to_numpy(cat_tautau), "costheta_tautau"] =  helicityCosTheta(ditau_cat_tautau, tau1_cat_tautau).to_numpy()
        
        cat_taue = ( (n_taus == 1) & (n_electrons == 1) )
        tau1_cat_taue = selectedTaus[cat_taue][:,0] 
        tau2_cat_taue = selectedElectrons[cat_taue][:,0] 
        ditau_cat_taue = tau1_cat_taue + tau2_cat_taue
        output.loc[ak.to_numpy(cat_taue), "costheta_tautau"] = helicityCosTheta(ditau_cat_taue, tau1_cat_taue).to_numpy()

        cat_taumu = ( (n_taus == 1) & (n_muons == 1) )
        tau1_cat_taumu = selectedTaus[(n_taus == 1) & (n_muons == 1)][:,0] 
        tau2_cat_taumu = selectedMuons[(n_taus == 1) & (n_muons == 1)][:,0] 
        ditau_cat_taumu = tau1_cat_taumu + tau2_cat_taumu
        output.loc[ak.to_numpy(cat_taumu), "costheta_tautau"] = helicityCosTheta(ditau_cat_taumu, tau1_cat_taumu).to_numpy()
        
        scale1fb = 1
        output["weight"] = 1
        if "Data" not in dataset:

            scale1fbfile = 'data/samples_and_scale1fb_HHggTauTau.json'
            import os.path
            if not os.path.isfile(scale1fbfile):
                scale1fbfile = scale1fbfile.replace("data", "/hadoop/cms/store/user/hmei/workflowtest/test_coffea_VBF")

            with open(scale1fbfile) as json_file:
                inputs = json.load(json_file)
    
            samplename = "_".join(dataset.split("_")[:-1])
            year = dataset.split("_")[-1]
            scale1fb = inputs[samplename][year]["metadata"]["scale1fb"]
            logger.debug(f"scale1fb for {samplename} is {scale1fb}")
            logger.debug(f'genWeight for {samplename} is {events["genWeight"][:10]}')
            output["weight"] = events["genWeight"]*scale1fb*self.lumis[year]
            logger.debug(f'event weight for {samplename} in year {year} is {events["genWeight"][0]*scale1fb*self.lumis[year]}')

        if "GJets" in dataset:
            events.metadata["dataset"] = events.metadata["dataset"].split("_")[0]
        output["dataset"] = events.metadata["dataset"]

        output["vbfjets_eta1_times_eta2"] = -999.0
        output["vbfjets_abs_eta_diff"] = -999.0
        output["vbfjets_mjj"] = -999.0
        two_selectedVBFJets = selectedVBFJets[n_vbf_jets >= 2]
        output.loc[ak.to_numpy(n_vbf_jets >= 2), "vbfjets_eta1_times_eta2"] = ak.to_numpy(two_selectedVBFJets[:, 0].eta * two_selectedVBFJets[:, 1].eta)
        output.loc[ak.to_numpy(n_vbf_jets >= 2), "vbfjets_abs_eta_diff"] = ak.to_numpy(abs( two_selectedVBFJets[:, 0].eta - two_selectedVBFJets[:, 1].eta) )
        output.loc[ak.to_numpy(n_vbf_jets >= 2), "vbfjets_mjj"] = ak.to_numpy( (two_selectedVBFJets[:, 0] + two_selectedVBFJets[:, 1]).mass )

        import time
        import random
        for key, task in get_worker().tasks.items():
            if task.state == "executing":
                job_id = f"{key[-32:]}_{time.time()}_{dataset}"

        localoutputname = f"out-{self.job_tag}-{job_id}.parquet"
        #output[ak.to_numpy(all_cuts)].to_parquet(f"./{localoutputname}")
        output[ak.to_numpy(all_cuts)].to_parquet(f"{self.output_path}/{localoutputname}")

        #cmd_gfal_copy = f"env -i X509_USER_PROXY=${{X509_USER_PROXY}} gfal-copy -p -f -t 4200 --verbose file://`pwd`/{localoutputname} gsiftp://gftp.t2.ucsd.edu/{self.output_path}/{localoutputname} --checksum ADLER32"
        #call(cmd_gfal_copy, shell=True)

        #logger.debug(f"output: {output.head()}")
        #logger.debug(f"mask: {all_cuts}, type of mask: {type(all_cuts)}")
        #return output[ak.to_numpy(all_cuts)]
        return 0
    
    def postprocess(self, accumulator):
        return super().postprocess(accumulator)
        #return accumulator
