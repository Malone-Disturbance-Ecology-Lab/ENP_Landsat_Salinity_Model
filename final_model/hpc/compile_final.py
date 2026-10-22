"""Compile the final one-step model on Bouchet (tl2cgen) and check it against scikit-learn. Usage: python compile_final.py <input_dir>"""
import sys, pickle, numpy as np, treelite, tl2cgen
IN = sys.argv[1]; M = pickle.load(open(f'{IN}/final_model_onestep.pkl', 'rb')); m = M['m3']
tl2cgen.export_lib(treelite.sklearn.import_model(m), toolchain='gcc', libpath=f'{IN}/m3_final.so', params={'parallel_comp': 8}, nthread=4)
X = np.random.default_rng(0).normal(size=(5000, len(M['F']))).astype('float32') * 50
d = float(np.max(np.abs(tl2cgen.Predictor(f'{IN}/m3_final.so', nthread=4).predict(tl2cgen.DMatrix(X)).ravel() - m.predict(X))))
print('max |compiled - sklearn|', d); assert d < 1e-4
