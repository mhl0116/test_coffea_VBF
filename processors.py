import awkward as ak
from coffea import hist, processor
import pandas as pd

# register our candidate behaviors
from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

class VBFHHggtautauProcessor(processor.ProcessorABC):
    '''
    https://github.com/CoffeaTeam/coffea/blob/65978ea299eee653f01ccfbc11af9f0339ffc7c2/tests/test_dask_pandas.py
    https://github.com/CoffeaTeam/coffea/blob/bbfd34414530f981529ccecf16e1585d293c5389/coffea/processor/test_items/NanoTestProcessorPandas.py
    '''
    def __init__(self):
        self._accumulator = processor.dict_accumulator({
            "sumw": processor.defaultdict_accumulator(float),
            "mass": hist.Hist(
                "Events",
                hist.Cat("dataset", "Dataset"),
                hist.Bin("mass", "$m_{\ta\ta}$ [GeV]", 60, 60, 120),
            ),
        })

    @property
    def accumulator(self):
        return self._accumulator

    def process(self, events):

        output = pd.DataFrame() 
        dataset = events.metadata['dataset'] # not clear why this is needed
        
        from utils import wrap_items
        item_names_to_be_wrapped = ["selectedPhoton", "tau", "electron", "muon", "jet"]
        photons, taus, electrons, muons, jets = wrap_items(events, item_names_to_be_wrapped)

        from utils import select_tau #, select_electron, select_muon, select_jet
        selectedPhotons = photons
        selectedMuons = select_muon(muons, selectedPhotons)
        selectedElectrons = select_electron(electrons, selectedPhotons)
        selectedTaus = select_tau(taus, selectedPhotons, selectedElectrons, selectedMuons)
        #selectedJets = select_jet(jets, selectedPhotons) 

        #print (f'photon: {photons[0][0]}')
        print (f'tau: {taus[0][0]}')
        print (f'tau: {selectedTaus}')
        #print (f'electron: {electrons[0][0]}')
        #print (f'muon: {muons[0][0]}')
        #print (f'jet: {jets[0][0]}')

        dipho = photons[:, 0] + photons[:, 1]

        #cut = (ak.num(taus) == 2) & (ak.sum(taus.charge) == 0)
        # add first and second tau in every event together
        #ditau = taus[cut][:, 0] + taus[cut][:, 1]

        #output["sumw"][dataset] += len(events)
        #output["mass"].fill(
        #    dataset=dataset,
        #    mass=dipho.mass,
        #)

        return output

    def postprocess(self, accumulator):
        return accumulator
