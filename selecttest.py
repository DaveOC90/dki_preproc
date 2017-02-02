import nipype.pipeline.engine as pe
from nipype.interfaces import afni,fsl,ants,dipy
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype
import os,glob,sys
from nipype.workflows.dmri.fsl.dti import create_eddy_correct_pipeline



## Intialize variables and workflow
globaldir='/home/davidoconner/dki_preproc/'
workdir='/home/davidoconner/dki_preproc/working_test'
preproc = pe.Workflow(name='preprocflow')
preproc.base_dir = workdir

ipdir='/home/davidoconner/hbnssi_rawdata/'
sublist=[g.split('/')[-2] for g in glob.glob(ipdir+'*/')]
seslist=['ses-SSV1']#list(set([g.split('/')[-2] for g in glob.glob(ipdir+'*/*/')]))

# Topup Acquisition Parameters
#acqparm='0 -1 0 0.0684\n0 1 0 0.0684'
#index=' '.join([1 for x in range(0,numvols)])

# BBR Schedule
bbrsched='/usr/share/fsl/5.0/etc/flirtsch/bbr.sch'
## Setup data managment nodes

infosource = pe.Node(util.IdentityInterface(fields=['subject_id','session_id']),name="infosource")
infosource.iterables = [('subject_id', sublist),('session_id', seslist)]

templates={
'anat' : ipdir+'{subject_id}/{session_id}/anat/{subject_id}_{session_id}_acq-MEMPRAGE_T1w.nii.gz'
}

selectfiles = pe.Node(nio.SelectFiles(templates,base_directory=ipdir),name="selectfiles")

datasink = pe.Node(nio.DataSink(base_directory=globaldir, container=workdir),name="datasink")



## Anat Preproc
# Skullstrip  MPRAGE
anat_skullstrip = pe.Node(interface=afni.preprocess.SkullStrip(),
                                  name='anat_skullstrip')
anat_skullstrip.inputs.args = '-o_ply'
anat_skullstrip.inputs.outputtype = 'NIFTI_GZ'


# Mask MPRAGE
anat_brain_only = pe.Node(interface=afni.preprocess.Calc(),
                        name='anat_brain_only')
anat_brain_only.inputs.expr = 'a*step(b)'
anat_brain_only.inputs.outputtype = 'NIFTI_GZ'



# FAST Node
segment=pe.Node(interface = fsl.FAST(), name='segment')
segment.inputs.img_type = 1
segment.inputs.segments = True
segment.inputs.probability_maps = True
segment.inputs.out_basename = 'segment'

# Split fast prob maps, picking wm
wmprob_split = pe.Node(interface=util.Select(), name = 'wmprob_split')
wmprob_split.inputs.index=[2]



wmmapbin = pe.Node(interface=fsl.UnaryMaths(),name='wmmapbin')
#fslmaths opname_fast_pve_2 -thr 0.5 -bin opname_fast_wmseg
# in_file -> out_file
wmmapbin.inputs.args='-thr 0.5 -bin'



preproc.connect([



    (infosource,selectfiles,[('subject_id', 'subject_id'),('session_id', 'session_id')]),

    (selectfiles, anat_skullstrip, [('anat','in_file')]),
    (selectfiles, anat_brain_only, [('anat','in_file_a')]),
    (anat_skullstrip, anat_brain_only, [('out_file', 'in_file_b')]),
                            
    (anat_brain_only, segment, [('out_file',  'in_files')]),
    (segment,wmprob_split, [('probability_maps','inlist')])#,
    #(wmprob_split,wmmapbin, (['out','in_file']))

])

preproc.run('MultiProc',plugin_args={'n_procs':4})
preproc.write_graph()