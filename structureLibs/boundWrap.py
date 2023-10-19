from orderParam_lib_dev import getBoundWrap

def func1(topFile, frame, 
          watInds, watHInds, 
          solInds, solHInds, solCInds, solOInds, solNInds, solSInds):
    boundInds, wrapInds, shellInds, nonShellInds = getBoundWrap(topFile, frame,
                                                                watInds, watHInds,
                                                                solInds, solHInds,
                                                                solCInds, solOInds,
                                                                solNInds, solSInds, 
                                                                cutoff = 4.6)

    subInds = [boundInds, wrapInds, shellInds, nonShellInds]
    return subInds
