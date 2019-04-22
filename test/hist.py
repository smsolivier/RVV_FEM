#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt

import subprocess 
import re 
from collections import defaultdict

exe = 'avl' 

subprocess.run(['make', '-j', exe])

binary = subprocess.run(['spike', '-g', '--isa=rv64gcv', 'pk', exe], 
	stderr=subprocess.PIPE)
hist = binary.stderr.splitlines()

dump = subprocess.run(['riscv64-unknown-elf-objdump', '-d', exe], 
	stdout=subprocess.PIPE) 
pcs = dump.stdout.splitlines() 

insns = {} 
for line in pcs: 
	line = line.decode('utf-8')
	try:
		insn = re.findall(r'\s\s\s[0-9a-fA-F]+:\t[0-9a-fA-F]+\s*(\w*)', line)[0]
		pc = re.findall(r'\s\s\s([0-9a-fA-F]+):', line)[0] 
		insns[pc] = insn 
	except:
		insn = '' 

count = defaultdict(int) 

for line in hist: 
	line = line.decode('utf-8')
	pc = line.split(' ')[0] 
	c = line.split(' ')[1] 
	if (pc in insns):
		count[insns[pc]] += int(c)
		
print(count['vadd'])
order = sorted(count.items(), key=lambda item: item[1], reverse=True)
print(order[:30])