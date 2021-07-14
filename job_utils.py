def make_fileset(samplelist, yearlist, use_xrootd):
    
    fileset = {}

    for sample in samplelist:
        for year in ["2016", "2017", "2018"]:
            if year not in yearlist: continue
            namelist = f"./data/{sample}_{year}.txt"
            fnamelist = open(namelist)
            fnames = fnamelist.read().splitlines()

            if use_xrootd:
                fnames = [fname.replace("/hadoop/cms","root://redirector.t2.ucsd.edu/") for fname in fnames]

            fileset[f"{sample}_{year}"] = fnames

    return fileset
