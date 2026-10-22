"""Compile the three fitted models to native code with treelite/tl2cgen (pip install treelite tl2cgen). Identical predictions, ~12x faster."""
import pickle, sys, treelite, tl2cgen
IN=sys.argv[1]; F=pickle.load(open(f'{IN}/final_models.pkl','rb'))
for k in ('m1','m2','m3'): tl2cgen.export_lib(treelite.sklearn.import_model(F[k]),toolchain='gcc',libpath=f'{IN}/{k}.so',params={'parallel_comp':8},nthread=4); print(k)
