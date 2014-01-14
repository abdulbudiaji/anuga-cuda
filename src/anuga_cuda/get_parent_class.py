import os
import inspect


def get_parent_class(classObj):
    assert inspect.isclass(classObj)
    return inspect.getmro(classObj)
                
