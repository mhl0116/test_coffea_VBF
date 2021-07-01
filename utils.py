import awkward as ak
# register our candidate behaviors
from coffea.nanoevents.methods import candidate
ak.behavior.update(candidate.behavior)

def wrap_items(events, item_names):

    wrapped_items = []
    if "selectedPhoton" in item_names:
        wrapped_items.append(
            ak.zip(
                {
                "pt": events["selectedPhoton_pt"],
                "eta": events["selectedPhoton_eta"],
                "phi": events["selectedPhoton_phi"],
                "mass": events["selectedPhoton_mass"],
                # photon has charge 0
                "charge": events["selectedPhoton_mass"] - events["selectedPhoton_mass"],
                },
                with_name="PtEtaPhiMCandidate",
                #behavior=candidate.behavior,
            )
        )

    if "tau" in item_names:
        wrapped_items.append(
            ak.zip(
                {
                "pt": events["Tau_pt"],
                "eta": events["Tau_eta"],
                "phi": events["Tau_phi"],
                "mass": events["Tau_mass"],
                "charge": events["Tau_charge"],
                "decay_mode": events["Tau_idDecayModeNewDMs"],
                "dz": events["Tau_dz"],
                "deepTau_vs_j": events["Tau_idDeepTau2017v2p1VSjet"],
                "deepTau_vs_e": events["Tau_idDeepTau2017v2p1VSe"],
                "deepTau_vs_m": events["Tau_idDeepTau2017v2p1VSmu"],
                },
                with_name="PtEtaPhiMCandidate",
                #behavior=candidate.behavior,
            )
        )

    if "electron" in item_names:
        wrapped_items.append(
            ak.zip(
                {
                "pt": events["Electron_pt"],
                "eta": events["Electron_eta"],
                "phi": events["Electron_phi"],
                "mass": events["Electron_mass"],
                "charge": events["Electron_charge"],
                "dxy": events["Electron_dxy"],
                "dz": events["Electron_dz"],
                "mvaId_wp90": events["Electron_mvaFall17V2Iso_WP90"],
                "mvaId_noIso_wp90": events["Electron_mvaFall17V2noIso_WP90"],
                "pfRelIso03": events["Electron_pfRelIso03_all"],
                },
                with_name="PtEtaPhiMCandidate",
                #behavior=candidate.behavior,
            )
        )

    if "muon" in item_names:
        wrapped_items.append(
            ak.zip(
                {
                "pt": events["Muon_pt"],
                "eta": events["Muon_eta"],
                "phi": events["Muon_phi"],
                "mass": events["Muon_mass"],
                "charge": events["Muon_charge"],
                "dxy": events["Muon_dxy"],
                "dz": events["Muon_dz"],
                #"mvaId": events["Muon_"],
                "pfRelIso03": events["Muon_pfRelIso03_all"],
                },
                with_name="PtEtaPhiMCandidate",
                #behavior=candidate.behavior,
            )
        )

    if "jet" in item_names:
        wrapped_items.append(
            ak.zip(
                {
                "pt": events["Jet_pt"],
                "eta": events["Jet_eta"],
                "phi": events["Jet_phi"],
                "mass": events["Jet_mass"],
                "charge": events["Jet_mass"] - events["Jet_mass"],
                "neEmEF": events["Jet_neEmEF"],
                "neHEF": events["Jet_neHEF"],
                "chHEF": events["Jet_chHEF"],
                "nConstituents": events["Jet_nConstituents"],
                "btagDeepFlavB": events["Jet_btagDeepFlavB"],
                "btagDeepFlavC": events["Jet_btagDeepFlavC"],
                },
                with_name="PtEtaPhiMCandidate",
                behavior=candidate.behavior,
            )
        )

    return wrapped_items


def select_tau(taus, selectedPhotons, selectedElectrons, selectedMuons):
    
    pt_cut = taus.pt > 18 
    eta_cut = abs(taus.eta) < 2.3 
    decay_mode_cut = taus.decay_mode == True
    dz_cut = abs(taus.dz) < 0.2 

    id_electron_cut = taus.deepTau_vs_e >= 1 # vvloose 
    id_muon_cut = taus.deepTau_vs_m >= 0 # vloose 
    id_jet_cut = taus.deepTau_vs_j >= 3 # vloose

    # if dR(tau, pho/ele/mu) < threshold, then taus.nearest(XXX, threshold) is its dR, otherwise None
    # ak.fill_none(, -1) makes the value for taus that will be saved -1 
    # if electron or muon is not present, dR_ele_cut is also -1, then this means the tau being considered will not be filtered, as there is nothing to calculated dR with
    dR_pho_cut =  ak.fill_none( taus.nearest( selectedPhotons, threshold = 0.2 ).pt, -1 ) < 0
    dR_ele_cut =  ak.fill_none( taus.nearest( selectedElectrons, threshold = 0.2 ).pt, -1 ) < 0
    dR_muon_cut =  ak.fill_none( taus.nearest( selectedMuons, threshold = 0.2 ).pt, -1 ) < 0

    #dR_pho_cut = object_selections.select_deltaR(events, taus, photons, options["taus"]["dR_pho"], debug)
    #dR_muon_cut = object_selections.select_deltaR(events, taus, muons, options["taus"]["dR_lep"], debug)
    #dR_ele_cut = object_selections.select_deltaR(events, taus, electrons, options["taus"]["dR_lep"], debug)

    tau_cut = pt_cut & eta_cut & decay_mode_cut & dz_cut & id_electron_cut & id_muon_cut & id_jet_cut & dR_pho_cut & dR_muon_cut & dR_ele_cut

    return taus[tau_cut]

def select_electron(electrons, selectedPhotons):
    
    pt_cut = electrons.pt > 7 
    eta_cut = abs(electrons.eta) < 2.5 
    dxy_cut = abs(electrons.dz) < 0.045 
    dz_cut = abs(electrons.dz) < 0.2 

    eleId_cut = (electrons.mvaId_wp90 == True | ((electrons.mvaId_noIso_wp90 == True) & (electrons.pfRelIso03 < 0.3)))
    dR_pho_cut =  ak.fill_none( electrons.nearest( selectedPhotons, threshold = 0.2 ).pt, -1 ) < 0
    mZ_cut = electrons 

    electron_cut = pt_cut & eta_cut & dxy_cut & dz_cut & eleId_cut & dR_pho_cut & mZ_cut

    return electrons[electron_cut]
