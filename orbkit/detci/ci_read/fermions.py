import numpy

from orbkit.display import display
from orbkit.qcinfo import CIinfo
from orbkit.read.tools import descriptor_from_file
from orbkit.units import ev_to_ha
from .tools import multiplicity

#################################################################################
#  TDSCF
#################################################################################
def fermions_tdscf(fname,td_typ='rpa',st_typ='singlet',select_state=None,threshold=0.0,**kwargs):
  '''Reads FermiONs++ TDDFT output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of TDDFT calculation from FermiONs++...')
  #
  # future user-input...
  #
  sst = 'Singlets'
  fspin = 1.0
  if st_typ == 'triplet':
    fspin = 3.0
    sst = 'Triplets'

  ss0 = ''
  if td_typ == 'tda':
    ss0 = 'Tamm-Dancoff Excitation Energies [%s]' % sst
  elif td_typ == 'rpa':
    ss0 = 'RPA Excitation Energies [%s]' % sst
  elif td_typ == 'stda':
    ss0 = 'Simplified Tamm-Dancoff Excitation Energies [%s]' % sst
  elif td_typ == 'srpa':
    ss0 = 'Simplified RPA Excitation Energies [%s]' % sst

  nocca = -1
  noccb = -1

 # Initialize variables
  ci = []
  ci_flag = False
  prttol = False
  init_state = False
  rhfspin = 0
  spin = 'Unknown'
  nel = 0
  deex = []
  
  if isinstance(select_state,int): select_state = [select_state]

  if isinstance(fname, str):
    filename = fname
    fname = descriptor_from_file(filename, index=0, ci_descriptor=True)
  else:
    filename = fname.name

  st = 0
  for line in fname:
    thisline = line.split()             # The current line split into segments
    #--- Check the file for keywords ---
    # Initialize Hartree-Fock ground state
    if 'No. of alpha-electrons:' in line:
      nocca = int(line.split()[3])
    elif 'No. of beta-electrons:' in line:
      noccb = int(line.split()[3])
    elif 'Final SCF energy:' in line:
        ci.append(CIinfo(method='tddft'))
        ci[-1].info   = []
        ci[-1].coeffs = []
        ci[-1].occ    = []
        ci[-1].occ.append([0,0])
        ci[-1].coeffs.append(1.0)
        ci[-1].info = {'state': '0',
                       'energy': float(thisline[3]),
		       'energy_nm': 0.0,
                       'fileinfo': filename,
                       'read_threshold': threshold,
                       'spin': spin,
		       'f_0i': 0.0}
    # Initialize new excited state
    if st == 0: # look for RPA
      if ss0 in line:
        st = 1
    elif st == 1:
      if 'Converged' in line and '==' in line:
        st = 2
    elif st == 2:
      if 'Excited State' in line: # and 'eV' in line and 'nm' in line:
        if select_state is None or int(thisline[2]) in select_state:
          init_state = True
          tddft_skip = 1
          ci.append(CIinfo(method='tddft'))
          ci[-1].info   = []
          ci[-1].coeffs = []
          ci[-1].occ    = []
          ci[-1].info = {'state': thisline[2],
                         'energy': 0.0, #float(thisline[-6])*ev_to_ha + ci[0].info['energy'],
          	       'energy_nm': 0.0, #float(thisline[-4]),
                         'fileinfo': filename,
                         'read_threshold': threshold,
                         'spin': fspin, #thisline[3].split('-')[0],
          	       'f_0i': 0.0} #float(thisline[8].replace('=',' ').split()[-1])}
          deex.append([])
      if init_state == True:
        if not tddft_skip:
          if 'Excitation Energy:' in line:
            ci[-1].info['energy'] = float(thisline[2])*ev_to_ha + ci[0].info['energy'] 
          elif ' nm' in line:
            ci[-1].info['energy_nm'] = float(thisline[0])
#                           'spin': thisline[3].split('-')[0],
          elif len(line.strip()) == 0: #thisline == [] or '->' not in line and '<-' not in line:
            init_state = False
          else:
            if '<-->' in line:
              if abs(float(thisline[-2])) > threshold:
                ex  = thisline[-2]
                dex = thisline[-2]
                thisline = line.split('<-->') #replace('->','-> ').split()
                s0 = thisline[0].split('(')[1].split(')')[0].strip()
                nocc = nocca
                if 'beta' in line:
                  nocc = noccb
                s1 = int(thisline[1].split('(')[1].split(')')[0].strip())+nocc
                tmp_occ = [s0,'%i'%s1]
                ci[-1].occ.append(tmp_occ)
                ci[-1].coeffs.append(float(ex)*numpy.sqrt(2))
                deex[-1].append(float(dex)*numpy.sqrt(2))
        elif tddft_skip:
          tddft_skip -= 1
      if '-------------------' in line or 'Excite Total:' in line:
        st = 3

  fname.close()
  deex = numpy.array(deex)
  #--- Calculating norm of CI states
  display('\nIn total, %d states have been read.' % len(ci)) 
  display('Norm of the states:')
  for i in range(len(ci)):
    j = numpy.array(ci[i].coeffs,dtype=float)
    norm = numpy.sum(j**2)
    ci[i].coeffs = j
    # Write Norm to log-file
    display('\tState %s:\tNorm = %0.8f (%d Coefficients)' % (ci[i].info['state'],norm, len(ci[i].coeffs)))
    # Transform to numpy arrays
    ci[i].occ = numpy.array([s for s in ci[i].occ],dtype=numpy.intc)-1
  
  return ci

def fermions_tddft(fname,select_state=None,threshold=0.0,**kwargs):
  '''Reads FermiONs++ TDDFT output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of TDDFT calculation from FermiONs++...')
  return fermions_tdscf(fname,td_typ='rpa',st_typ='singlet',select_state=select_state,threshold=threshold,**kwargs)

def fermions_tda_tddft(fname,select_state=None,threshold=0.0,**kwargs):
  '''Reads FermiONs++ TDDFT output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of TDA-TDDFT calculation from FermiONs++...')
  return fermions_tdscf(fname,td_typ='tda',st_typ='singlet',select_state=select_state,threshold=threshold,**kwargs)

def fermions_srpa(fname,select_state=None,threshold=0.0,**kwargs):
  '''Reads FermiONs++ TDDFT output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of simplified RPA calculation from FermiONs++...')
  return fermions_tdscf(fname,td_typ='srpa',st_typ='singlet',select_state=select_state,threshold=threshold,**kwargs)

def fermions_stda(fname,select_state=None,threshold=0.0,**kwargs):
  '''Reads FermiONs++ TDDFT output. 
  
  **Parameters:**
  
    fname: str, file descriptor
      Specifies the filename for the input file.
      fname can also be used with a file descriptor instad of a filename.
    select_state : None or list of int, optional
      If not None, specifies the states to be read (0 corresponds to the ground 
      state), else read all electronic states.
    threshold : float, optional
      Specifies a read threshold for the CI coefficients.
  
  **Returns:**
  
    ci : list of CIinfo class instances
      See :ref:`Central Variables` for details.
  '''
  display('\nReading data of simplified TDA calculation from FermiONs++...')
  return fermions_tdscf(fname,td_typ='stda',st_typ='singlet',select_state=select_state,threshold=threshold,**kwargs)

#################################################################################
#  CASCI -> todo...
#################################################################################
