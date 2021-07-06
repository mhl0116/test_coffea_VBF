import uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema
#uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource
### https://uproot.readthedocs.io/en/latest/uproot.source.xrootd.html
uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.XRootDSource

import coffea.processor as processor
from processors import VBFHHggtautauProcessor

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')

file_handler = logging.FileHandler('logs/example.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)


samplelist = ["HH_ggTauTau"] #, "VBF_HH_ggTauTau"]
import job_utils
fileset = job_utils.make_fileset(samplelist, use_xrootd=True) 
logger.debug(f"fileset: {fileset}")

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
        client = Client(memory_limit='2GB', n_workers=8, threads_per_worker=1)

        out = processor.run_uproot_job(
            fileset,
            treename = 'Events',
            processor_instance = VBFHHggtautauProcessor(),
            #executor=processor.futures_executor,
            #executor_args={"schema": None, "workers": 3, "use_dataframes": True}, # our skim only works with None if we want to use selectedPhoton..
            executor=processor.dask_executor,
            executor_args={"schema": None, "client": client, "use_dataframes": True},
            )

        logger.debug(f"columns: {out.columns}")
        logger.debug(f"head of df: {out.head()}")

        import dask.dataframe as dd
        dd.to_parquet(out, path="./outputs/test2.parquet")

        client.shutdown()

if __name__ == '__main__':
    run(False)
