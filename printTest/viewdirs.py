import numpy as np
import nibabel as nb
import dipy.reconst.dti as dti
from dipy.segment.mask import median_otsu
from dipy.reconst.dti import fractional_anisotropy, color_fa, lower_triangular
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table
from dipy.data import get_sphere
from dipy.viz import fvtk

def gradientpicture(bvalpath, bvecpath, niipath):

	img=nb.Nifti1Image.load(niipath)
	bvals,bvecs=read_bvals_bvecs(bvalpath,bvecpath)
	gtab = gradient_table(bvals, bvecs)
	data = img.get_data()
	print('data.shape (%d, %d, %d, %d)' % data.shape)

	maskdata, mask = median_otsu(data, 3, 1, True,  dilate =2) #vol_idx=range(10, 50), dilate=2)
	print('maskdata.shape (%d, %d, %d, %d)' % maskdata.shape)

	print('Creating and fitting tensor model')
	tenmodel = dti.TensorModel(gtab)

	tenfit = tenmodel.fit(maskdata)
	print('Computing anisotropy measures (FA, MD, RGB)')
	print('FA first')	
	FA = fractional_anisotropy(tenfit.evals)

	FA[np.isnan(FA)] = 0

	fa_img = nb.Nifti1Image(FA.astype(np.float32), img.get_affine())
	nb.save(fa_img, 'tensor_fa.nii.gz')

	evecs_img = nb.Nifti1Image(tenfit.evecs.astype(np.float32), img.get_affine())
	nb.save(evecs_img, 'tensor_evecs.nii.gz')

	print('Now MD')
	MD1 = dti.mean_diffusivity(tenfit.evals)
	nb.save(nb.Nifti1Image(MD1.astype(np.float32), img.get_affine()), 'tensors_md.nii.gz')
	MD2 = tenfit.md

	print('Now RGB')
	FA = np.clip(FA, 0, 1)
	RGB = color_fa(FA, tenfit.evecs)
	nb.save(nb.Nifti1Image(np.array(255 * RGB, 'uint8'), img.get_affine()), 'tensor_rgb.nii.gz')

	print('Computing tensor ellipsoids in a part of the splenium of the CC')

	
	sphere = get_sphere('symmetric724')

	
	ren = fvtk.ren()

	#evals = tenfit.evals[13:43, 44:74, 28:29]
	#evecs = tenfit.evecs[13:43, 44:74, 28:29]

	indstr=":,60:61,:"

	evals = eval("tenfit.evals["+indstr+"]")
	evecs = eval("tenfit.evecs["+indstr+"]")

	evals=np.reshape(evals,(evals.shape[0],evals.shape[2],1,3))
	evecs=np.reshape(evecs,(evecs.shape[0],evecs.shape[2],1,3,3))


	#cfa = RGB[13:43, 44:74, 28:29]
	cfa = eval("RGB["+indstr+"]")

	cfa=np.reshape(cfa,(cfa.shape[0],cfa.shape[2],1,3))

	cfa /= cfa.max()

	fvtk.add(ren, fvtk.tensor(evals, evecs, cfa, sphere))

	print('Saving illustration as tensor_ellipsoids.png')
	fvtk.record(ren, n_frames=1, out_path='tensor_ellipsoids.png', size=(1200, 1200))

	fvtk.clear(ren)

	#tensor_odfs = tenmodel.fit(data[20:50, 55:85, 38:39]).odf(sphere)
	#tensor_odfs = tenmodel.fit(data[20:50, 55:85, 38:39]).odf(sphere)

	data = eval("data["+indstr+"]")

	tensor_odfs = tenmodel.fit(data).odf(sphere)
 	tensor_odfs = tenmodel.fit(data).odf(sphere)
	fvtk.add(ren, fvtk.sphere_funcs(tensor_odfs, sphere, colormap=None))
	print('Saving illustration as tensor_odfs.png')
	fvtk.record(ren, n_frames=1, out_path='tensor_odfs.png', size=(600, 600))


	## Look at 3dautobox 
if __name__ == '__main__':
	
	import sys

	niipath=sys.argv[1]
	bvalpath=sys.argv[2]
	bvecpath=sys.argv[3]

	gradientpicture(bvalpath, bvecpath, niipath)