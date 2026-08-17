import subprocess

prog_loc1 = 'gc_content/'

cmd1 = (
		f'python3 {prog_loc1}GC_content.py ../data/smallgenes/ --apc_results '
		f'../data/APCisos_base/ ../data/APCisos_optiso/ '
		f'../data/APCisos_nmd/ ../data/APCisos_optiso_nmd/'
	).split(' ')

cmd2 = (
		f'python3 make_all_svg_files.py ../data/smallgenes/ '
		f'../data/APCisos_base/ --out_dir apc_base_gcc_svgs/ '
		f'--gcc gcc_res/rnaseq_gcc/ gcc_res/APC_base_gcc/'
	).split(' ')

cmd3 = (
		f'python3 make_all_svg_files.py ../data/smallgenes/ '
		f'../data/APCisos_optiso/ --out_dir apc_optiso_gcc_svgs/ '
		f'--gcc gcc_res/rnaseq_gcc/ gcc_res/APC_optiso_gcc/'
	).split(' ')
	
cmd4 = (
		f'python3 make_all_svg_files.py ../data/smallgenes/ '
		f'../data/APCisos_nmd/ --out_dir apc_nmd_gcc_svgs/ '
		f'--gcc gcc_res/rnaseq_gcc/ gcc_res/APC_nmd_gcc/'
	).split(' ')
	
cmd5 = (
		f'python3 make_all_svg_files.py ../data/smallgenes/ '
		f'../data/APCisos_optiso_nmd/ --out_dir '
		f'apc_optiso_nmd_gcc_svgs/ --gcc gcc_res/rnaseq_gcc/ '
		f'gcc_res/APC_optiso_nmd_gcc/'
	).split(' ')
	
print('calculate gc content')
subprocess.run(cmd1)
print('make base svgs')
subprocess.run(cmd2)
print('make optiso svgs')
subprocess.run(cmd3)
print('make nmd svgs')
subprocess.run(cmd4)
print('make optiso_nmd svgs')
subprocess.run(cmd5)
