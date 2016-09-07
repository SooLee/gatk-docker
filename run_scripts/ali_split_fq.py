#!/usr/bin/python
# -*- coding: utf-8 -*-
#### ali_split_fq.py
#### made by Min-Seok Kwon
#### 2016-09-06 14:21:13
#########################
import sys
import os
from datetime import datetime
import argparse

VERSION = "v0.01"
VERSION_DATE = "08/15/2016"
PROG = sys.argv[0].split("/")[-1].strip()
PROG_PATH = sys.argv[0]
#TITLE = os.getcwd().split("/")[-1]
OPT = ""
FQLINE = 20*4000000
MAXMEM_MCORE = 4
BINSIZE = 10*1000*1000

def fileSave (path, cont, opt, gzip_flag = "n"):
	if gzip_flag == "gz":
		import gzip
		f = gzip.open(path, opt)
		f.write(cont)
		f.close()
	else:
		f = file (path, opt)
		f.write(cont)
		f.close

def comma(value):
	return "{:,}".format(value)

def nf(k,lenk = 4):
	sk = str(k)
	for i in range(len(sk), lenk):
		sk = "0" + sk
	return sk

def step1_split_fq():
	global OPT

	now1 = datetime.now()
	now3 = now1

	j = 0
	j2 = 0  ### each turn
	j3 = 0	### each lane
	k = 0
	lane_k = {}
	lane_readno = {}
	fmap = {}

	pre_lane = ""
	logfile = OPT.wd + "/" + OPT.readRno + ".fq.log"
	
	for line in sys.stdin:
		if line[0] == "@" and j % 4 == 0:
			if len(line.split(":")) > 3:
				lane = "_".join(line[1:].split(":")[:4])
			else:
				lane = line[1:].split(".")[0] ### @SRR3184307.1 1/1
			#print lane, pre_lane



			try:
				tmp = fmap[lane]
				if lane_readno[lane] % (FQLINE) == 0:
					fmap[lane].close()
					lane_k[lane] += 1
					fqname = OPT.wd + "/" +lane+"."+nf(lane_k[lane],3)+"."+ OPT.readRno+".fq"
					fmap[lane] = open(fqname,"w")	

					i = j / 4
					i2 = j2 / 4
					i3 = j3 / 4
					now2 = datetime.now()
					t1 = now2-now1
					t2 = now2-now3
					log=[str(datetime.now()),"turn "+comma(i2)+ " reads", "lane " +comma(lane_readno[lane])+ " reads","total " +comma(i)+ " reads", str(t2.total_seconds()) + " sec/turn", str(t1.total_seconds()) + " sec"] 
					print "\t".join(log)
					now3 = now2				

				lane_readno[lane] += 1
			except KeyError:
				lane_readno[lane] = 1
				try:
					k = lane_k[lane]
				except KeyError:
					lane_k[lane] = 1

				fqname = OPT.wd + "/" +lane+"."+nf(lane_k[lane],3)+"."+ OPT.readRno+".fq"
				fmap[lane] = open(fqname,"w")


		
		fmap[lane].write(line)
		j += 1
		j2 += 1
		j3 += 1

	for lane in fmap.keys():
		fmap[lane].close()


def dispatch_job():
	step1_split_fq()

def init():
	global OPT

	p1 = argparse.ArgumentParser(usage='%(prog)s <sub-command> [options]', description='%(prog)s '+VERSION+"("+VERSION_DATE+")"+': alignment pipeline manager')
	p1.add_argument('-readRno', dest='readRno', default="R1", choices=['R1','R2'], help='read paid ID')	
	p1.add_argument('-wd', dest='wd', default=False, help="working directory")

	
	#p1.add_argument('-step9', dest='step9', action="store_true", default=False, help="mkdup")
	#p1.add_argument('-fq', dest='fq', default="", help="fq")
	#p1.add_argument('-key', dest='key', default="", help="key")
	#p1.add_argument('-mergelane', dest='mergelane', action="store_true", default=False, help="merge lane")

	if len(sys.argv) == 1 or (len(sys.argv) == 2 and sys.argv[1][0] != '-'):
		sys.argv.append('-h')
	OPT = p1.parse_args()
	dispatch_job()	

if __name__ == "__main__":
	init()
