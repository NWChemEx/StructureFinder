import os
import parallelzone as pz
import sys
import unittest



if __name__ == '__main__':
    rv = pz.runtime.RuntimeView()

    my_dir = os.path.dirname(os.path.realpath(__file__))
    root_dir = os.path.dirname(os.path.dirname(os.path.dirname(my_dir)))
    src_dir  = os.path.join(root_dir, 'src', 'python')
    sys.path.append(src_dir)

    loader = unittest.TestLoader()
    tests  = loader.discover(my_dir)
    testrunner = unittest.runner.TextTestRunner()
    ret = not testrunner.run(tests).wasSuccessful()
    sys.exit(ret)

