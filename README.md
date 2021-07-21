To set up environment

```
conda create uproot dask dask-jobqueue matplotlib pandas jupyter hdfs3 pyarrow fastparquet numba numexpr --prefix ./envs
conda install -c conda-forge jupyterlab

conda activate /afs/cern.ch/work/h/hmei/myWorkspace/envs
conda config --set env_prompt '(Hgg_{name})'
conda install -c conda-forge boost-histogram
conda run pip install jupyter-server-proxy coffea autopep8 jupyter_nbextensions_configurator klepto yahist

conda activate /afs/cern.ch/work/h/hmei/myWorkspace/envs
conda pack --prefix ~/myWorkspace/Hgg/envs -f --format tar.gz --compress-level 9 -j 8 --exclude "*pyc" --exclude "*.js.map" --exclude "*a" --exclude "*pandoc"
mv envs.tag.gz make_dask_cluster/Hgg_envs.tar.gz
```
Then put the  file in the make\_dask\_cluster directory
Scripts from make\_dask\_cluster foler are copied and modified from this [repo](https://github.com/aminnj/daskucsd).

I usually create a cluster in a seperate notebook or ipython, as in createcluster.ipynb, then copy the address of client to example.py:

```
client = Client("tcp://169.228.130.37:12318")
```
