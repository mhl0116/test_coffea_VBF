import json

import uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema
uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource

import coffea.processor as processor
from processors import VBFHHggtautauProcessor

with open('data/samples_and_scale1fb_HHggTauTau.json') as json_file:
    inputs = json.load(json_file)

fnamelist = open("./data/VBF_HH_ggTauTau_2018.txt")
fnames = fnamelist.read().splitlines()

use_xrootd = True
if use_xrootd:
    fnames = [fname.replace("/hadoop/cms","root://redirector.t2.ucsd.edu/") for fname in fnames]

fileset = {"VBF_HH_ggTauTau": fnames[:5]}

def run(useNanoEvents):

    if useNanoEvents == True:

        file = uproot.open(fnames[1])
        events = NanoEventsFactory.from_root(
            file,
            entry_stop=10000,
            metadata={"dataset": "VBFHHggtautau"},
            schemaclass=NanoAODSchema,
        ).events()
        p = MyProcessor()
        out = p.process(events)
        return out, events

    else:

        from dask.distributed import Client
        client = Client(memory_limit='2GB', n_workers=3, threads_per_worker=1)

        #p = MyProcessor()
        out = processor.run_uproot_job(
            fileset,
            treename = 'Events',
            processor_instance = VBFHHggtautauProcessor(),
            executor=processor.futures_executor,
            executor_args={"schema": None, "workers": 3}, # our skim only works with None if we want to use selectedPhoton..
            #executor=processor.dask_executor,
            #executor_args={"schema": None, "client": client},
            )
        return out

if __name__ == '__main__':
    run(False)
