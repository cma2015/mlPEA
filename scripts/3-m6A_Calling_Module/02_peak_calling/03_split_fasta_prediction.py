import numpy as np
from ushuffle import shuffle
import sys
import re

# 获取输入文件路径
input_fa = sys.argv[1]
output_path = sys.argv[2]

with open(f'{input_fa}', 'r') as f1:
    pos_lines = f1.readlines()

for i in range(0, len(pos_lines)):
    pos_lines[i] = pos_lines[i].rstrip('\n')
    pos_lines[i] = pos_lines[i].lstrip('>')

def  count_list(std:list, tongji):
    from collections import Counter
    name = Counter()
    for num in std:
        name[num] += 1
    print(name)
    
dic = {
  'A': 0, 'a': 0,
  'C': 1, 'c': 1,
  'G': 2, 'g': 2,
  'T': 3, 't': 3,
}

def asc2one(a):
    a = np.copy(a)  # 将数组转换为可写
    length = len(a)
    for i in range(0, length):
        a[i] = dic[chr(a[i])]
    return a
  
token = []
labels = []

for i in range(1, len(pos_lines), 2):
    cleaned_line = re.sub("N", "", pos_lines[i])
    sub_array = asc2one(np.frombuffer(cleaned_line.encode(), dtype=np.int8))
    token.append(sub_array)
    labels.append([1])

token = np.array(token, dtype=object)
labels = np.array(labels, dtype=np.int8)

np.save(f'{output_path}/test_token.npy', token)
np.save(f'{output_path}/test_label.npy', labels)