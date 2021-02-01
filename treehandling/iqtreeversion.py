import os, re
def find_iqtree_version():
	versionpattern = r'version (\w+)'
	command = os.popen("iqtree -h").read()
	command = command.split("\n")[0]
	version = re.search(versionpattern, command)
	if version:
		iqversion = version.group(1)
		print(command, "; thus version", iqversion)
	else:
		iqversion = 1
	return iqversion

version = find_iqtree_version()
