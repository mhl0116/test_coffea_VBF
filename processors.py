import awkward as ak
from coffea import hist, processor
import pandas as pd
import json

# register our candidate behaviors
from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')

file_handler = logging.FileHandler('logs/processor.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

class VBFHHggtautauProcessor(processor.ProcessorABC):
    '''
    https://github.com/CoffeaTeam/coffea/blob/65978ea299eee653f01ccfbc11af9f0339ffc7c2/tests/test_dask_pandas.py
    https://github.com/CoffeaTeam/coffea/blob/bbfd34414530f981529ccecf16e1585d293c5389/coffea/processor/test_items/NanoTestProcessorPandas.py
    '''
    lumis = {"2016": 35.9, "2017": 41.9, "2018": 59.7}

    def __init__(self):
        pass

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):

        output = pd.DataFrame() 
        dataset = events.metadata['dataset'] # identifier of data/MC 
        
        from utils import wrap_items
        item_names_to_be_wrapped = ["selectedPhoton", "tau", "electron", "muon", "jet"]
        photons, taus, electrons, muons, jets = wrap_items(events, item_names_to_be_wrapped)

        ### object selection
        from utils import select_tau, select_electron, select_muon, select_jet
        selectedPhotons = photons
        selectedMuons = select_muon(muons, selectedPhotons)
        selectedElectrons = select_electron(electrons, selectedPhotons)
        selectedTaus = select_tau(taus, selectedPhotons, selectedElectrons, selectedMuons)
        selectedJets = select_jet(jets, selectedPhotons, selectedMuons, selectedElectrons, selectedTaus) 
        events["selectedPhoton"] = selectedPhotons
        events["selectedMuon"] = selectedMuons
        events["selectedElectron"] = selectedElectrons
        events["selectedTau"] = selectedTaus
        events["selectedJet"] = selectedJets

        dipho = selectedPhotons[:, 0] + selectedPhotons[:, 1]
        events["diphoton"] = dipho 

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
        all_cuts = os_cut | one_tau_zero_lep_cut

        #### save relevant info
        output["run"] = ak.to_numpy(events["run"])
        output["event"] = ak.to_numpy(events["event"])
        output["dataset"] = events.metadata["dataset"]

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
        
        scale1fb = 1
        output["weight"] = 1
        if "Data" not in dataset:
            with open('data/samples_and_scale1fb_HHggTauTau.json') as json_file:
                inputs = json.load(json_file)
            samplename = "_".join(dataset.split("_")[:-1])
            year = dataset.split("_")[-1]
            scale1fb = inputs[samplename][year]["metadata"]["scale1fb"]
            logger.debug(f"scale1fb for {samplename} is {scale1fb}")
            logger.debug(f'genWeight for {samplename} is {events["genWeight"][:10]}')
            output["weight"] = events["genWeight"]*scale1fb*self.lumis[year]
            logger.debug(f'event weight for {samplename} in year {year} is {events["genWeight"][0]*scale1fb*self.lumis[year]}')

        #logger.debug(f"output: {output.head()}")
        #logger.debug(f"mask: {all_cuts}, type of mask: {type(all_cuts)}")
        return output[ak.to_numpy(all_cuts)]

    def postprocess(self, accumulator):
        return accumulator
