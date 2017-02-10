import os, h5py

def get_ccs_passes(filename):
    hf = h5py.File(filename)
    basename = os.path.basename(filename)[:-9]

    suffix = filename[-9:]
    assert suffix in ('.1.ccs.h5', '.2.ccs.h5', '.3.ccs.h5')

    for hole_num, num_passes in zip(hf['/PulseData/ConsensusBaseCalls/ZMW/HoleNumber'].value, hf['/PulseData/ConsensusBaseCalls/Passes/NumPasses'].value):
        if num_passes > 0:
            print '%s/%d/ccs' % (basename, hole_num), num_passes

if __name__ == '__main__':
    import sys
    for filename in sys.argv[1:]:
        get_ccs_passes(filename)