import os
import re

root_func = "evolve"

def get_class_method_hierarchy(file_name=None, file_string=None, file_line_list=None):
    code_line_lsit = []
    func_dic = {}
    if file_name is not None:
        if os.path.exists(file_name):
            with open(file_name, "r") as f:
                code_line_lsit = f.readlines()
    elif file_string:
        code_line_lsit = file_string.split("\n")
    elif file_line_list:
        code_line_lsit = file_line_list

    find_func(code_line_lsit, func_dic)
    for func in func_dic:
        find_subfunc(func, code_line_lsit, func_dic)

    
    root_list = [func for func in func_dic if len(func_dic[func][3]) and len(func_dic[func][4]) == 0]
    print root_list
    leaf_list = [func for func in func_dic if len(func_dic[func][3]) == 0]
    print leaf_list
    return leaf_list
    #for root in root_list:
    #    build_func_tree(root, func_dic)

    #for func in func_dic:
    #    print func, func_dic[func][2:]



def find_func(code_list, func_dic):
    pre_func = ""
    for i in range(len(code_list)):
        line = code_list[i]
        if re.findall("def.*\(", line):
            func_name = re.split('[\t ]*\([\t ]*', re.split('[\t ]*def[\t ]*', line)[1])[0]
            # Code beginning line #, ending #, height, subfunction list, parentfunction list
            func_dic[func_name] = [i, -1, -1, [], []]
            if func_dic.has_key(pre_func):
                func_dic[pre_func][1] = i
            pre_func = func_name

    func_dic[pre_func][1] = len(code_list) - 1



def find_subfunc(parent_function, code_list, func_dic):
    line_no_B, line_no_E = func_dic[parent_function][0:2]
    subfunc_list = func_dic[parent_function][3]
    for line in code_list[line_no_B: line_no_E]:
        if re.findall("\w* *\(", line):
            subfunc_name_list = [re.split(" *\(", tmp_subfunc_name)[0] for tmp_subfunc_name in re.findall("[a-zA-Z_]\w* *\(", line)]
            for subfunc_name in subfunc_name_list:
                if func_dic.has_key(subfunc_name) and subfunc_name != parent_function:
                    subfunc_list.append(subfunc_name)
                    if parent_function not in func_dic[subfunc_name][4]:
                        func_dic[subfunc_name][4].append(parent_function)
                    
            if len(func_dic[parent_function][3]):
                func_dic[parent_function][2] = 0



def build_func_tree(parent_func, func_dic):
    func_height = 0
    for subfunc in func_dic[parent_func][3]:
        tmp_func_height = build_func_tree(subfunc, func_dic)
        if func_height <= tmp_func_height:
            func_height = tmp_func_height + 1
    func_dic[parent_func][2] = func_height
    return func_height



tab_space = 4

def indentation_analysier(code_line, level, calling_func):
    tab_num = 0 
    global tab_space
    for i in range(level * tab_space):
        if code_list[i] == '\t':
            tab_num += 1
        elif code_line[i: i + tab_space] == ' '*tab_space:
            tab_num += 1
        else: 
            break
    return level - tab_num 



if __name__ == "__main__":
    get_class_method_hierarchy("gpu_domain_advanced.py")
