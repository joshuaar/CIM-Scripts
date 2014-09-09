import numpy as np

"""
Simulates a basic model of peptide arrays based on first order kinetics
"""

def showPareto(m, a):
	s = -np.random.pareto(a, 1000) + m 
	s = s - min(s) + 10**-7
	import matplotlib.pyplot as plt
	count, bins, ignored = plt.hist(s, 100)
	fit = a*m**a/bins**(a+1)
	#plt.plot(bins, max(count)*fit/max(fit),linewidth=2, color='r')
	print min(s), max(s)
	plt.show()

class KineticProps(object):
	"""
	Intrinsic and kinetic properties of peptide array system
	params:
		kOn: On rate, scaler, s^-1 units
		kOff: Off rate matrix, matrix, nAbs X nPeptides
		kI: Irreversible rate, scaler, optional (defaults to 0) 
		cDist: Concentration distribution, nAbs x 1
	"""
	def __init__(self,kOn,kOff,cDistAb,cDistPep,kI = 0):
		self.kOn = kOn #scaler
		self.kOff = np.array(kOff) #matrix
		self.kI = kI #scaler
		self.cDistAb = np.array(cDistAb)
		self.cDistPep = np.array(cDistPep)

	@classmethod
	def fromPareto(cls, nAbs, nPeps, minOff = 10**-7, a = 485, kOn = 0.1 kI = 0.01):
		"Determines off rates for each peptide/antibody interaction using Pareto distribution"
		samps = -np.random.pareto(a, nPeps*nAbs)
		samps = samps - min(samps) + minOff #adjust to proper kOff units
		kOff = np.reshape(samps,(nAbs,nPeps))
		cDistAb = [1./nAbs for i in range(nAbs)] #assume concentrations are uniformly distributed
		cdistPep = [1./nPeps for i in range(nPeps)] #same with peptide "concentrations"
		return cls(kOn,kOff,cDistAb,kI)

	def mod_bulkPeptide(self, kOff = 10**-2, prop=2.0/3):
		"""
		Adds a "bulk peptide" to represent the interstitial space.
		Creates a new column in the kOff matrix, and adjusts the cDistPep distribution
		params:
			kOff: the off rate for the new interstitial peptide
			prop: the proportion of total peptides that are interstitial peptide
		"""
		kOffCol = np.array([[kOff for i in range(self.kOff.shape[0])]]).transpose()
		self.kOff = np.hstack((self.kOff,kOffCol)) #modify the kOff matrix
		self.cDistPep = np.append(self.cDistPep * (1.0 - prop), prop) #modify the distribution of peptide conc
		#now we've included an interstitial peptide as last entry in kOff and cDistPep

class AssayProps(object):
	"""
	Assay Properties of peptide array system
	params:
		cAb: total concentration of antibody in nM (scaler)
		cPep: total concentration of peptide in nM (scaler)
	"""
	def __init__(self,cAb,cPep):
		self.cAb = cAb
		self.cPep = cPep
	@classmethod
	def fromVols(cls,vol,d,cAb,surface=0.49)
		"""
		Calculate Assay properties from densities and volumes
		params:
			vol: total volume in uL
			d: peptide density in nm^-2 (n peptides per square nanometer)
			cAb: the total antibody concentration
			surface: surface area in cm^2
		"""
		dcm = d * 10**14 #peptide density in square centimeters
		npep = dcm * surface #number of peptides
		A = 6.02 * 10**23 #avagadro's number
		molesPep = npep/A
		volInL = (vol * 1.0) / 10**6
		cPep = molesPep / volInL #concentration of peptide
		return cls(cAb=cAb,cPep=cPep)

class SimulateSystem(object):
	"""
	Simulates system of KineticProps and AssayProps
	params:
		KineticProps
		AssayProps
	"""
	def __init__(self, KineticProps, AssayProps):
		self.AssayProps = AssayProps
		self.KineticProps = KineticProps
		self.cAbsFree, self.cPepsFree, self.cAbsBound, self.cPepsBound, self.cAbsBoundIrrev = self.calcInitial(AssayProps,KineticProps)

	def simulate(self, t=1, dt = 0.01)
		"Simulates t hours using dt seconds as step"
		stateMatrix = []
		nSteps = t * 60 * 60 * dt
		for i in range(nSteps):
			stateMatrix.append(self.__simulateOne(dt))
		return np.array(stateMatrix)
			
	
	def __simulateOne(self, dt):
		kI = self.KineticProps.kI
		kOn = self.KineticProps.kOn
		kOff = self.KineticProps.kOff #off rate matrix
		for i, abRow in enumerate(kOff): #enumerate over antibodies in kOff
			for j, pepCol in enumerate(abRow):
				cAbsBoundT_1 = self.cAbsBound[i][j]
				cAbsFreeT_1 = self.cAbsFree[i]
				self.cAbsFree[i] += (kOff[i][j]*cAbsBoundT_1 - kOn*cAbsFreeT_1*self.cPepsFree[j])*dt
				self.cAbsBound[i][j] += (kOn*cAbsfreeT_1*self.cPepsFree[j] - (kI + kOff[i][j])*self.cAbsBound[i][j])*dt 
				self.cAbsBoundIrrev[i][j] += (kI*cAbsBoundT_1)*dt
				self.cPepsFree[j] += (kOff[i][j]*cAbsBoundT_1 - kOn*cAbsfreeT_1*self.cPepsFree[j])*dt
				#This simulates one step of the model
		return self.getSystemState()

	def getSystemState(self):
		"Get the state of the system in some datastructure"
		pep_flourescence = np.zeros(self.cPepsFree.shape)
		for i, abRow in enumerate(self.cAbsBound):
			for j,pepCol in enumerate(abRow):
				pep_flourescence += self.cAbsBound[i][j]
				pep_flourescence += self.cAbsBoundIrrev[i][j]
		return pep_flourescence

	def calcInitial(self,AssayProps,KineticProps):
		"Calculate the initial conditions"
		cAb = AssayProps.cAb
		cPep = AssayProps.cPep
		cDistAb = KineticProps.cDistAb
		cDistPep = KineticProps.cDistPep
		cAbsFree = cAb * cDistAb
		cPepsFree = cPep * cDistPep
		cAbsBound = np.zeros(KineticProps.kOff.shape)
		cAbsBoundIrrev = np.zeros(KineticProps.kOff.shape)
		cPepsBound = np.zeros(len(cPepsFree))
		return cAbsFree, cPepsFree, cAbsBound, cPepsBound, cAbsBoundIrrev

if __name__ == "__main__":
	"CLI (not finished)"
        parser = argparse.ArgumentParser()
        parser.add_argument("infile")
        parser.add_argument("-c", "--coeff", action="store_true")
        parser.add_argument("-a", "--alphabet", action="store_true")
        parser.add_argument("-i", "--intercept", action="store_true", default=False)
        parser.add_argument("-m", "--remove_end", action="store_true")
        parser.add_argument("-s", "--shuffle", action="store_true")
        parser.add_argument("-v", "--variance_inflation", action="store_true")
        parser.add_argument("-o", "--output")
        parser.add_argument("-t", "--log_transform", action="store_true")
        parser.add_argument("-l","--limit", nargs=2, type=float)
        parser.add_argument("-e","--exclude", help="exclude amino acids")
        parser.add_argument("--col", type=int)
        args = parser.parse_args()

