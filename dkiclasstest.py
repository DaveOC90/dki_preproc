
from nipype.interfaces import dipy
dkifit = dipy.DKI()
dkifit.inputs.in_file='/home/davidoconner/dki_preproc/working/preprocflow/_session_id_ses-SSV10_subject_id_sub-0031121/eddycorrect/eddy_corrected.nii.gz'
dkifit.inputs.in_bval='/home/davidoconner/hbnssi_rawdata/sub-0031121/ses-SSV1/dwi/sub-0031121_ses-SSV1_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.bval'
dkifit.inputs.in_bvec='/home/davidoconner/hbnssi_rawdata/sub-0031121/ses-SSV1/dwi/sub-0031121_ses-SSV1_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.bvec'
dkifit.run()
