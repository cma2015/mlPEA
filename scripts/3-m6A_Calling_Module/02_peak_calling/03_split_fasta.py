import numpy as np
from ushuffle import shuffle
import sys
import re
from sklearn.model_selection import train_test_split

with open(sys.argv[1] + '/pos.fa', 'r') as f1:
    pos_lines = f1.readlines()

for i in range(0, len(pos_lines)):
    pos_lines[i] = pos_lines[i].rstrip('\n')
    pos_lines[i] = pos_lines[i].lstrip('>')


with open(sys.argv[1] + '/neg.fa', 'r') as f2:
    neg_lines = f2.readlines()


for i in range(0, len(neg_lines)):
    neg_lines[i] = neg_lines[i].rstrip('\n')
    neg_lines[i] = neg_lines[i].lstrip('>')

def count_list(std:list,tongji):
    from collections import Counter
    name=Counter()
    for  num in std:
        name[num]+=1
    print(name)
    
dic = {
  'A' : 0, 'a' : 0,
  'C' : 1, 'c' : 1,
  'G' : 2, 'g' : 2,
  'T' : 3, 't' : 3,
}

def asc2one(a):
    length = len(a)
    for i in range(0, length):
        a[i] = dic[chr(a[i])]
    return a
  
token = []
lables = []

# for i in range(1, len(pos_lines), 2):
#     sub_array = asc2one(np.frombuffer(re.sub("N", "", pos_lines[i]), dtype=np.int8))
#     token.append(sub_array)
#     lables.append([1])
# 
# for i in range(1, len(neg_lines), 2):
#     sub_array = asc2one(np.frombuffer(re.sub("N", "", neg_lines[i]), dtype=np.int8))
#     token.append(sub_array)
#     lables.append([0])

for i in range(1, len(pos_lines), 2):
    cleaned_line = re.sub("N", "", pos_lines[i]).encode('utf-8')
    sub_array = np.frombuffer(cleaned_line, dtype=np.int8).copy()  # 复制到一个可写数组
    sub_array = asc2one(sub_array)
    token.append(sub_array)
    lables.append([1])

for i in range(1, len(neg_lines), 2):
    cleaned_line = re.sub("N", "", neg_lines[i]).encode('utf-8')
    sub_array = np.frombuffer(cleaned_line, dtype=np.int8).copy()  # 复制到一个可写数组
    sub_array = asc2one(sub_array)
    token.append(sub_array)
    lables.append([0])

# 将 token 转换为 NumPy 数组，并指定数据类型为 object，以便处理不同长度的子数组
all_tokens = np.array(token, dtype=object)
all_labels = np.array(lables, dtype=np.int8)

train_token,test_token,train_label,test_label = train_test_split(all_tokens,all_labels,test_size=0.2,random_state=0)
test_token,valid_token,test_label,valid_label = train_test_split(test_token,test_label,test_size=0.5,random_state=0)

print(train_token.shape,test_token.shape,valid_token.shape,train_label.shape,test_label.shape,valid_label.shape)

np.save(sys.argv[1] + '/data/train_token.npy', train_token)
np.save(sys.argv[1] + '/data/test_token.npy', test_token)
np.save(sys.argv[1] + '/data/valid_token.npy', valid_token)
np.save(sys.argv[1] + '/data/train_label.npy', train_label)
np.save(sys.argv[1] + '/data/test_label.npy', test_label)
np.save(sys.argv[1] + '/data/valid_label.npy', valid_label)