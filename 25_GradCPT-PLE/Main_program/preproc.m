%% Philippdata preproc


Indir = '/BICNAS2/group-northoff/yujia/gradcptdata/BIDS';
Outdir = '/BICNAS2/group-northoff/yujia/gradcptdata/BIDS_derivative';

sub_dir = dir(fullfile(Indir,'0*'));

sub_label = [];
sub_num = length(sub_dir);

for i = 1:10
    sub_label = sprintf('%s %s',sub_label,sub_dir(i).name);
end

CommandInit = sprintf('singularity run --cleanenv -B /BICNAS2:/BICNAS2 /home/test-oc/singularity/images/fmriprep-23.0.2.simg %s %s participant --fs-license-file %s',Indir,Outdir,'/home/yao/Downloads/matlab_toolbox/Freesurfer_license/license.txt');

Command = sprintf('%s --skip-bids-validation',CommandInit);
Command = sprintf('%s --participant-label %s',Command,sub_label);
Command = sprintf('%s --output-spaces fsaverage5',Command);
Command = sprintf('%s --fs-subjects-dir %s',Command,fullfile(Outdir,'freesurfer'));
Command = sprintf('%s --cifti-output 91k',Command);
Command = sprintf('%s --nthreads 16',Command);
Command = sprintf('%s --low-mem',Command);

system(Command,'-echo')

