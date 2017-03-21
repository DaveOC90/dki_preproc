import nipype.pipeline.engine as pe
from nipype.interfaces import afni,fsl,ants,dipy
import nipype.interfaces.io as nio
import nipype.interfaces.utility as util
import nipype
import os,glob,sys
from nipype.workflows.dmri.fsl.dti import create_eddy_correct_pipeline




def hex2float(inmat):

    '''
    Function to convert matrices in hexadecimal format to floating point format
    Input: double space delimited .mat file from FSL
    Output: double space delimited .mat file
    '''

    # Read in data as list of lists, seperating strings based on double spaces and newline characters
    data=[l.strip().split('  ') for l in open(inmat,'rU')]
    # convert each element to float
    datanew=[map(lambda x: float.fromhex(x),d) for d in data]
    # turn list of lists back into one long string
    opdatanew='\n'.join(['  '.join(map(str,d)) for d in datanew])
    # write mat to same directory as input mat
    opname=inmat.replace('.mat','_float.mat')
    # write string
    fo=open(opname,'w')
    fo.write(opdatanew)
    fo.close()

    return opname






def create_anat_preproc(name='anat_preproc'):
    ## Anat Preproc
    '''
    Setting up workflow for skullstripping, segmenting, and spatially
    normalizing a T1 weighted MR Image
    '''

    anat_preproc = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['anat']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['ss_brain',
                                                        'wm_map_bin',
                                                        'mni_xfm']),
                         name='outputspec')


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
    wmprob_split.inputs.index=2

    # Binarizing white matter probability map
    wmmapbin = pe.Node(interface=fsl.maths.MathsCommand(),name='wmmapbin')
    wmmapbin.inputs.args='-thr 0.5 -bin'


    # Calculate ANTs Warp
    calculate_ants_warp = pe.Node(interface=ants.Registration(),
                name='calculate_ants_warp')
    calculate_ants_warp.inputs. \
        fixed_image = '/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain.nii.gz'
    calculate_ants_warp.inputs. \
        dimension = 3
    calculate_ants_warp.inputs. \
        use_histogram_matching=[ True, True, True ]
    calculate_ants_warp.inputs. \
        winsorize_lower_quantile = 0.01
    calculate_ants_warp.inputs. \
        winsorize_upper_quantile = 0.99
    calculate_ants_warp.inputs. \
        metric = ['MI','MI','CC']
    calculate_ants_warp.inputs. \
        metric_weight = [1,1,1]
    calculate_ants_warp.inputs. \
        radius_or_number_of_bins = [32,32,4]
    calculate_ants_warp.inputs. \
        sampling_strategy = ['Regular','Regular',None]
    calculate_ants_warp.inputs. \
        sampling_percentage = [0.25,0.25,None]
    calculate_ants_warp.inputs. \
        number_of_iterations = [[1000,500,250,100], \
        [1000,500,250,100], [100,100,70,20]]
    calculate_ants_warp.inputs. \
        convergence_threshold = [1e-8,1e-8,1e-9]
    calculate_ants_warp.inputs. \
        convergence_window_size = [10,10,15]
    calculate_ants_warp.inputs. \
        transforms = ['Rigid','Affine','SyN']
    calculate_ants_warp.inputs. \
        transform_parameters = [[0.1],[0.1],[0.1,3,0]]
    calculate_ants_warp.inputs. \
        shrink_factors = [[8,4,2,1],[8,4,2,1],[6,4,2,1]]
    calculate_ants_warp.inputs. \
         smoothing_sigmas = [[3,2,1,0],[3,2,1,0],[3,2,1,0]]
    calculate_ants_warp.inputs. \
         sigma_units = ['vox','vox','vox']
    calculate_ants_warp.inputs. \
        output_warped_image = True
    calculate_ants_warp.inputs. \
        output_inverse_warped_image = True
    calculate_ants_warp.inputs. \
        output_transform_prefix = 'xfm'
    calculate_ants_warp.inputs. \
        write_composite_transform = True
    calculate_ants_warp.inputs. \
    collapse_output_transforms = False


    anat_preproc.connect([

    (inputspec, anat_skullstrip, [('anat','in_file')]),
    (inputspec, anat_brain_only, [('anat','in_file_a')]),
    (anat_skullstrip, anat_brain_only, [('out_file', 'in_file_b')]),
                               
    (anat_brain_only, calculate_ants_warp, [('out_file',  'moving_image')]),
    (anat_brain_only, segment, [('out_file',  'in_files')]),
    (segment,wmprob_split, [('probability_maps','inlist')]),
    (wmprob_split,wmmapbin, [('out','in_file')]),

    (anat_brain_only, outputspec,[('out_file', 'ss_brain')]),
    (wmmapbin, outputspec, [('out_file', 'wm_map_bin')]),
    (calculate_ants_warp, outputspec, [('composite_transform', 'mni_xfm')])

    ])

    return anat_preproc


def create_diff_preproc(name='diff_preproc'):

    ## Diffusion Tensor Computation
    hexmat2fltmat = util.Function(input_names=["inmat"], output_names=["opmat"],function=hex2float)


    diff_preproc = pe.Workflow(name=name)
    
    inputspec = pe.Node(util.IdentityInterface(fields=['diff',
                                                       'bvals',
                                                       'bvecs',
                                                       'b0ap',
                                                       'b0pa']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['diff_preproc']),
                         name='outputspec')


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
    # fslmaths my_hifi_b0 -Tmean my_hifi_b0
    #output is 'file'

    bet_b0 = pe.Node(interface=fsl.BET(), name='bet_b0')
    bet_b0.inputs.mask = True
    #bet my_hifi_b0 my_hifi_b0_brain -m



    eddycorrect = pe.Node(interface=fsl.Eddy(), name='eddycorrect')
    eddycorrect.inputs.in_index = '/home/davidoconner/git/dki_preproc/index.txt'
    eddycorrect.inputs.in_acqp  = '/home/davidoconner/git/dki_preproc/acqparams.txt'
    eddycorrect.threads = 2
    #eddy --imain=data --mask=my_hifi_b0_brain_mask --acqp=acqparafslviewms.txt --index=index.txt --bvecs=bvecs --bvals=bvals --topup=my_topup_results --out=eddy_corrected_data


    ## Convert eddy mat to float
    h2f = pe.Node(interface=hexmat2fltmat,name='h2f')





    diff_preproc.connect([

    (inputspec,fslroi,[('diff','in_file')]),
    (fslroi, bet, [('roi_file', 'in_file')]),

    (inputspec,mergelistb0,[('b0ap','in1')]),
    (inputspec,mergelistb0,[('b0pa','in2')]),
    

    (mergelistb0,b0merge,[('out','in_files')]),
    (b0merge,topup,[('merged_file','in_file')]),

    (topup,b0_corrected_mean,[('out_corrected','in_file')]),
    (b0_corrected_mean,bet_b0,[('out_file','in_file')]),

    (inputspec, eddycorrect, [('bvals', 'in_bval')]),
    (inputspec, eddycorrect, [('bvecs', 'in_bvec')]),
    (inputspec, eddycorrect, [('diff','in_file')]),
    (bet_b0, eddycorrect, [('mask_file', 'in_mask')]),
    (topup,eddycorrect, [('out_fieldcoef','in_topup_fieldcoef')]),
    (topup,eddycorrect, [('out_movpar','in_topup_movpar')]),
    (eddycorrect,outputspec,[('out_corrected','diff_preproc')])

    ])

    return diff_preproc


def create_diff_norm(name='diff_norm'):

   
    diff_norm = pe.Workflow(name=name)
    
    hexmat2fltmat = util.Function(input_names=["inmat"], output_names=["opmat"],function=hex2float)

    inputspec = pe.Node(util.IdentityInterface(fields=['diff_preproc',
                                                       'ss_brain',
                                                       'wm_map_bin',
                                                       'mni_xfm']),
                        name='inputspec')

    outputspec = pe.Node(util.IdentityInterface(fields=['diff_norm']),
                         name='outputspec')

    ## COnvert b0 flirt to anat initialization mat to float
    h2f2 = pe.Node(interface=hexmat2fltmat,name='h2f2')


    fslroi_b0corr = pe.Node(interface=fsl.ExtractROI(), name='fslroi_b0corr')
    fslroi_b0corr.inputs.t_min = 0
    fslroi_b0corr.inputs.t_size = 1

    ## B0 to Anat Initial mat
    linear_reg_b0_init = pe.Node(interface=fsl.FLIRT(), name='linear_reg_b0_init')
    linear_reg_b0_init.inputs.cost = 'corratio'
    linear_reg_b0_init.inputs.dof = 6
    linear_reg_b0_init.inputs.interp = 'trilinear'


    ## B0 to Anat
    linear_reg_b0 = pe.Node(interface=fsl.FLIRT(), name='linear_reg_b0')
    linear_reg_b0.inputs.cost = 'bbr'
    linear_reg_b0.inputs.dof = 6
    linear_reg_b0.inputs.interp = 'nearestneighbour'


    ## Apply XFM
    app_xfm_lin = pe.Node(interface=fsl.ApplyXfm(),
                         name='app_xfm_lin')
    app_xfm_lin.inputs.apply_xfm = True

    
    ## Func to MNI
    b0_t1_to_mni = pe.Node(interface=ants.ApplyTransforms(), name='b0_t1_to_mni')
    b0_t1_to_mni.inputs.dimension = 3
    b0_t1_to_mni.inputs.reference_image='/usr/share/fsl/5.0/data/standard/MNI152_T1_3mm_brain.nii.gz'
    b0_t1_to_mni.inputs.invert_transform_flags = [False]
    b0_t1_to_mni.inputs.interpolation = 'NearestNeighbor'
    b0_t1_to_mni.inputs.input_image_type = 3


    diff_norm.connect([

    (eddycorrect,fslroi_b0corr,[('out_corrected','in_file')]),

    (fslroi_b0corr, linear_reg_b0_init, [('roi_file', 'in_file')]),
    (anat_brain_only, linear_reg_b0_init,[('out_file', 'reference')]),

    (linear_reg_b0_init, h2f2, [('out_matrix_file','inmat')]),

    (fslroi_b0corr, linear_reg_b0, [('roi_file', 'in_file')]),
    (anat_brain_only, linear_reg_b0,[('out_file', 'reference')]),
    (h2f2, linear_reg_b0,[('opmat', 'in_matrix_file')]),
    (wmmapbin, linear_reg_b0, [('out_file', 'wm_seg')]),

    (eddycorrect, app_xfm_lin, [('out_corrected','in_file')]),
    (anat_brain_only, app_xfm_lin, [('out_file','reference')]),

    (linear_reg_b0, h2f, [('out_matrix_file','inmat')]),
    (h2f, app_xfm_lin, [('opmat','in_matrix_file')]),

    (calculate_ants_warp, b0_t1_to_mni, [('composite_transform', 'transforms')]),
    (app_xfm_lin, b0_t1_to_mni, [('out_file','input_image')]),
    (b0_t1_to_mni,outputspec, [('output_image','diff_norm')]),

    ])

    return diff_norm

def create_tensor_model(name='tensor_model'):
    dtifit = pe.Node(interface=fsl.DTIFit(), name='dtifit')

def create_kurtosis_model(name='kurtosis_model'):
    dkifit = pe.Node(interface=dipy.DKI(), name='dkifit')



if __name__ == '__main__':

    ## Intialize variables and workflow
    globaldir='/home/davidoconner/dki_preproc/subflows/'
    workdir='/home/davidoconner/dki_preproc/subflows/working/'
    preproc = pe.Workflow(name='preprocflow')
    preproc.base_dir = workdir

    ipdir='/home/davidoconner/hbnssi_rawdata/'
    sublist=[g.split('/')[-2] for g in glob.glob(ipdir+'*/')]
    seslist=['ses-SSV1']#list(set([g.split('/')[-2] for g in glob.glob(ipdir+'*/*/')]))

    # Topup Acquisition Parameters
    # acqparm='0 -1 0 0.0684\n0 1 0 0.0684'
    #index=' '.join([1 for x in range(0,numvols)])

    # BBR Schedule
    bbrsched='/usr/share/fsl/5.0/etc/flirtsch/bbr.sch'
    ## Setup data managment nodes

    infosource = pe.Node(util.IdentityInterface(fields=['subject_id','session_id']),name="infosource")
    infosource.iterables = [('subject_id', sublist),('session_id', seslist)]

    templates={
    'dki' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.nii.gz', \
    'bvals' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.bval', \
    'bvecs' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DKI64DIRECTIONSAP3WEIGHTSAX_dwi.bvec', \
    'b0ap' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DWIB0APAX_dwi.nii.gz', \
    'b0pa' : ipdir+'{subject_id}/{session_id}/dwi/{subject_id}_{session_id}_acq-DWIB0PAAX_dwi.nii.gz', \
    'anat' : ipdir+'{subject_id}/{session_id}/anat/{subject_id}_{session_id}_acq-MEMPRAGE_T1w.nii.gz'
    }

    selectfiles = pe.Node(nio.SelectFiles(templates,base_directory=ipdir),name="selectfiles")

    datasink = pe.Node(nio.DataSink(base_directory=globaldir, container=workdir),name="datasink")

    diff_preproc=create_diff_preproc()


    preproc.connect([

        (infosource,selectfiles,[('subject_id', 'subject_id'),('session_id', 'session_id')]),
        
        (selectfiles,diff_preproc,[('dki','inputspec.diff')]),
        (selectfiles,diff_preproc,[('b0ap','inputspec.b0ap')]),
        (selectfiles,diff_preproc,[('b0pa','inputspec.b0pa')]),
        
        (selectfiles, diff_preproc, [('bvals', 'inputspec.bvals')]),
        (selectfiles, diff_preproc, [('bvecs', 'inputspec.bvecs')])

        #(eddycorrect, dtifit, [('out_corrected', 'dwi')]),
        #(infosource, dtifit, [('subject_id', 'base_name')]),
        #(bet, dtifit, [('mask_file', 'mask')]),
        #(selectfiles, dtifit, [('bvals', 'bvals')]),
        #(selectfiles, dtifit, [('bvecs', 'bvecs')]),

        #(dtifit,datasink,[('FA','@.dtifit.FA')]),
        #(dtifit,datasink,[('L1','@.dtifit.L1')]),
        #(dtifit,datasink,[('L2','@.dtifit.L2')]),
        #(dtifit,datasink,[('L3','@.dtifit.L3')]),
        #(dtifit,datasink,[('MD','@.dtifit.MD')]),
        #(dtifit,datasink,[('MO','@.dtifit.MO')]),
        #(dtifit,datasink,[('S0','@.dtifit.S0')]),
        #(dtifit,datasink,[('V1','@.dtifit.V1')]),
        #(dtifit,datasink,[('V2','@.dtifit.V2')]),
        #(dtifit,datasink,[('V3','@.dtifit.V3')]),
        #(dtifit,datasink,[('tensor','@.dtifit.tensor')]),




        #(eddycorrect, dkifit, [('out_corrected', 'in_file')]),
        #(selectfiles, dkifit, [('bvals', 'in_bval')]),
        #(selectfiles, dkifit, [('bvecs', 'in_bvec')]),

        #(dkifit,datasink,[('fa','@.dkimodel.FA')]),
        #(dkifit,datasink,[('md','@.dkimodel.MD')]),
        #(dkifit,datasink,[('rd','@.dkimodel.RD')]),
        #(dkifit,datasink,[('ad','@.dkimodel.AD')]),
        #(dkifit,datasink,[('mk','@.dkimodel.MK')]),
        #(dkifit,datasink,[('ak','@.dkimodel.AK')]),
        #(dkifit,datasink,[('rk','@.dkimodel.RK')])


    ])


    preproc.run('MultiProc',plugin_args={'n_procs':4})
    preproc.write_graph()