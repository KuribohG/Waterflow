import os
import sys
os.system("chmod +w /unsullied/sharefs/wangmengdi/wangmengdi_EMC")

if __name__ == "__main__":
	boxname = sys.argv[1]
	basename = boxname.split("/")[-1].rstrip(".box")
	dumppref = "objs/%s"%basename
	os.system("mkdir -p %s"%dumppref)
	print boxname,dumppref
	os.system("rlaunch --cpu=1 --gpu=0 --memory=32768 -- ./FluidSimulation %s %s"%(boxname,dumppref))
