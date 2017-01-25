import nipype.pipeline.engine as pe
from nipype.interfaces import afni,fsl
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype
import os,glob,sys
from nipype.workflows.dmri.fsl.dti import create_eddy_correct_pipeline

## Intialize variables and workflow
globaldir='/home/davidoconner/dki_preproc/'
workdir='/home/davidoconner/dki_preproc/working'
preproc = pe.Workflow(name='preprocflow')
preproc.base_dir = workdir

ipdir='/home/davidoconner/hbnssi_rawdata/'
sublist=[g.split('/')[-2] for g in glob.glob(ipdir+'*/')]
seslist=list(set([g.split('/')[-2] for g in glob.glob(ipdir+'*/*/')]))

# Topup Acquisition Parameters
#acqparm='0 -1 0 0.0684\n0 1 0 0.0684'
#index=' '.join([1 for x in range(0,numvols)])


## Setup data managment nodes

infosource = pe.Node(util.IdentityInterface(fields=['subject_id','session_id']),name="infosource")
infosource.iterables = [('subject_id', sublist),('session_id', seslist)]

templates={
'dki' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.nii.gz', \
'bvals' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.bval', \
'bvecs' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.bvec', \
'b0ap' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DWIB0APAX_dwi.nii.gz', \
'b0pa' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DWIB0PAAX_dwi.nii.gz'
}

selectfiles = pe.Node(nio.SelectFiles(templates,base_directory=ipdir),name="selectfiles")

datasink = pe.Node(nio.DataSink(base_directory=globaldir, container=workdir),name="datasink")

## Diffusion Tensor Computation


fslroi = pe.Node(interface=fsl.ExtractROI(), name='fslroi')
fslroi.inputs.t_min = 0
fslroi.inputs.t_size = 1

bet = pe.Node(interface=fsl.BET(), name='bet')
bet.inputs.mask = True
bet.inputs.frac = 0.34


mergelistb0 = pe.Node(interface=util.Merge(2), name='mergelistb0')

b0merge = pe.Node(interface=fsl.Merge(), name='b0merge')
b0merge.inputs.dimension='t'
# input 'in_files' list of items
# output merged_file

topup = pe.Node(interface=fsl.TOPUP(), name='topup')
topup.inputs.config='b02b0.cnf'
topup.inputs.encoding_file='/home/davidoconner/git/dki_preproc/acqparams.txt'


b0_corrected_mean = pe.Node(interface=fsl.maths.MeanImage(), name='b0_corrected_mean')
b0_corrected_mean.inputs.dimension='T'
#fslmaths my_hifi_b0 -Tmean my_hifi_b0
#output is 'file'

bet_b0 = pe.Node(interface=fsl.BET(), name='bet_b0')
bet_b0.inputs.mask = True
#bet my_hifi_b0 my_hifi_b0_brain -m



eddycorrect = pe.Node(interface=fsl.Eddy(), name='eddycorrect')
eddycorrect.inputs.in_index = '/home/davidoconner/git/dki_preproc/index.txt'
eddycorrect.inputs.in_acqp  = '/home/davidoconner/git/dki_preproc/acqparams.txt'
eddycorrect.threads = 2
#eddy --imain=data --mask=my_hifi_b0_brain_mask --acqp=acqparams.txt --index=index.txt --bvecs=bvecs --bvals=bvals --topup=my_topup_results --out=eddy_corrected_data


dtifit = pe.Node(interface=fsl.DTIFit(), name='dtifit')

preproc.connect([

	(infosource,selectfiles,[('subject_id', 'subject_id'),('session_id', 'session_id')]),
    (selectfiles,fslroi,[('dki','in_file')]),
    (fslroi, bet, [('roi_file', 'in_file')]),
    (selectfiles,mergelistb0,[('b0ap','in1')]),
    (selectfiles,mergelistb0,[('b0pa','in2')]),
    (mergelistb0,b0merge,[('out','in_files')]),
    (b0merge,topup,[('merged_file','in_file')]),

    (topup,b0_corrected_mean,[('out_corrected','in_file')]),
    (b0_corrected_mean,bet_b0,[('out_file','in_file')]),

    (selectfiles, eddycorrect, [('bvals', 'in_bval')]),
    (selectfiles, eddycorrect, [('bvecs', 'in_bvec')]),
    (selectfiles, eddycorrect, [('dki','in_file')]),
    (bet_b0, eddycorrect, [('mask_file', 'in_mask')]),
    (topup,eddycorrect, [('out_fieldcoef','in_topup_fieldcoef')]),
    (topup,eddycorrect, [('out_movpar','in_topup_movpar')]),

    (eddycorrect, dtifit, [('out_corrected', 'dwi')]),
    (infosource, dtifit, [('subject_id', 'base_name')]),
    (bet, dtifit, [('mask_file', 'mask')]),
    (selectfiles, dtifit, [('bvals', 'bvals')]),
    (selectfiles, dtifit, [('bvecs', 'bvecs')]),

    (dtifit,datasink,[('FA','@dtifitFA')]),
    (dtifit,datasink,[('L1','@dtifitL1')]),
    (dtifit,datasink,[('L2','@dtifitL2')]),
    (dtifit,datasink,[('L3','@dtifitL3')]),
    (dtifit,datasink,[('MD','@dtifitMD')]),
    (dtifit,datasink,[('MO','@dtifitMO')]),
    (dtifit,datasink,[('S0','@dtifitS0')]),
    (dtifit,datasink,[('V1','@dtifitV1')]),
    (dtifit,datasink,[('V2','@dtifitV2')]),
    (dtifit,datasink,[('V3','@dtifitV3')]),
    (dtifit,datasink,[('tensor','@dtifittensor')])
])

preproc.run('MultiProc',plugin_args={'n_procs':4})
preproc.write_graph()