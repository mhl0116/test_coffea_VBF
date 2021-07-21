import uproot
from coffea.nanoevents import NanoEventsFactory, BaseSchema, NanoAODSchema
#uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.MultithreadedXRootDSource
### https://uproot.readthedocs.io/en/latest/uproot.source.xrootd.html
uproot.open.defaults["xrootd_handler"] = uproot.source.xrootd.XRootDSource

import coffea.processor as processor
from processors import VBFHHggtautauProcessor
import job_utils

import logging

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')

file_handler = logging.FileHandler('logs/example.log')
file_handler.setFormatter(formatter)

logger.addHandler(file_handler)

def run(fileset, outpath, jobtag, client, useNanoEvents):

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


        out = processor.run_uproot_job(
            fileset,
            treename = 'Events',
            processor_instance = VBFHHggtautauProcessor(outpath, jobtag),
            #executor=processor.futures_executor,
            #executor_args={"schema": None, "workers": 10}, # our skim only works with None if we want to use selectedPhoton..
            executor=processor.dask_executor,
            executor_args={"schema": None, "client": client, "use_dataframes": True},
            )

        #logger.debug(f"columns: {out.columns}")
        #logger.debug(f"head of df: {out.head()}")

        #import dask.dataframe as dd
        #dd.to_parquet(out, path="./outputs/test2_useNumbaForDr.parquet")
        #dd.to_parquet(out, path="./outputs/test2.parquet")
        #dd.to_parquet(out, path="./outputs/" + outdfname)

        #client.shutdown()

if __name__ == '__main__':

    #samplelist = ["HH_ggTauTau", "VBF_HH_ggTauTau", "Data"]
    #samplelist = ["HH_ggTauTau"]
    samplelist = ["HH_ggTauTau", "DiPhoton", "GJets_HT40To100", "GJets_HT-100To200", "GJets_HT-200To400", "GJets_HT-400To600", "GJets_HT-600ToInf", "ZGamma", "VH"]
    fileset = job_utils.make_fileset(samplelist, ["2016","2017","2018"], use_xrootd=True) 

    #outdfname = "test_VBFyields_withdata.parquet" 
    outpath = "/hadoop/cms/store/user/hmei/workflowtest/dask_coffea"
    jobtag = "test_heliticy"

    from dask.distributed import Client
    client = Client(memory_limit='2GB', n_workers=10, threads_per_worker=1)
    #client = Client("tcp://127.0.0.1:30055")

    localfiles = ["./processors.py", "./utils.py"]

    for localfile in localfiles:
        client.upload_file(localfile)

    logger.debug(f"fileset: {fileset}")

    run(fileset=fileset, outpath=outpath, jobtag=jobtag, client=client, useNanoEvents=False)

    #client.close()
