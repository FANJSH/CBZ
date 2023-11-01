# GNU C++
CXX = g++ 
CXXFLAGS = -O4 -Wall
#CXXFLAGS = -Wall

# Intel C++
#CXX = icpc
#CXXFLAGS = -O3

AR = ar -crsv
RM = rm -f

TARGET = libCBG.a
OBJS = 	Numeric.o FunctionTable.o GeneralOption.o GroupData.o CartCore.o CartMeshInfo.o \
       	GeneralMesh.o Medium.o Bessel.o GeomVector.o GeomPlane.o IrregularGeometryInformation.o \
       	SphericalHarmonics.o GroupDataSet.o Nuclide.o \
       	DelayedNeutronData.o GeneralSystem.o \
       	UNC_Covariance.o \
       	UNC_CrossSection.o UNC.o UNC_Parameters.o UNC_Sensitivity.o \
       	CartesianLattice.o \
       	LibData.o GammaLibData.o GammaTool.o \
       	PLOS_system.o  \
       	PLOSE_system.o \
       	SNT_system.o Quadrature.o \
       	SNR_system.o SNR_quadrature.o  \
       	SNRZ_system.o SNRZ_quadrature.o \
       	PJI_trajectory.o PJI_trajectoryset.o PJI_system.o PJI_slabpij.o \
       	MEC_system.o \
       	PNX_system.o \
       	ABEMIE_core.o ABEMIE_medium.o ABEMIE_plane.o ABEMIE_region.o \
       	DHEX_system.o \
       	OnePointCalculator.o SelfShieldingCalculator.o \
       	Burnup.o BurnupChain.o Burner.o FRBurnerRZ.o Cooler.o \
       	CrossSectionEdit.o \
       	HigherModeCalculator.o \
       	FRDesignTool.o FRDesign.o \
       	PKR.o \
       	OktavianTool.o SVDTool.o ENDFIDTranslator.o MATIDTranslator.o \
       	BurnupChainGenerator.o \
       	KineticUtility.o \
       	PincellCalculator.o \
       	RandomSamplingSupporter.o \
       	SensitivityData.o SUTool.o \
       	YieldDecayCovariance.o \
       	MultiStepCalculation.o \
       	MulticellBurner.o \
       	GeneralBurner.o \
       	Pij_dyn.o \
       	FMG_Tool.o \
		GData.o \
		HeatConductionCalculator.o jsmest.o jsmest_wrapper.o TwoPhaseFlowCalculator.o \
		RMforFR.o BNCTPerturbator.o\
		SomeOutputMethod.o\
		ResonanceCalculationFan.o

OUTPUTS =
BACKUPS = *~

all: ${TARGET}

${TARGET}: ${OBJS}
	${AR} ${TARGET} ${OBJS}

## -- for COMMON Library --

Numeric.o:Numeric.cxx Numeric.h 
	${CXX} ${CXXFLAGS} -c Numeric.cxx

FunctionTable.o:FunctionTable.cxx FunctionTable.h Numeric.h
	${CXX} ${CXXFLAGS} -c FunctionTable.cxx

GeneralOption.o:GeneralOption.cxx GeneralOption.h Numeric.h
	${CXX} ${CXXFLAGS} -c GeneralOption.cxx

GroupData.o:GroupData.cxx GroupData.h Numeric.h MultiStepCalculation.h
	${CXX} ${CXXFLAGS} -c GroupData.cxx

MultiStepCalculation.o:MultiStepCalculation.cxx MultiStepCalculation.h GroupData.h Numeric.h
	${CXX} ${CXXFLAGS} -c MultiStepCalculation.cxx

CartCore.o:CartCore.cxx CartCore.h Numeric.h
	${CXX} ${CXXFLAGS} -c CartCore.cxx

CartMeshInfo.o:CartMeshInfo.cxx CartMeshInfo.h CartCore.h Numeric.h
	${CXX} ${CXXFLAGS} -c CartMeshInfo.cxx

GeneralMesh.o:GeneralMesh.cxx GeneralMesh.h Medium.h Numeric.h
	${CXX} ${CXXFLAGS} -c GeneralMesh.cxx

Medium.o:Medium.cxx Medium.h GroupData.h GroupDataSet.h Nuclide.h Numeric.h MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c Medium.cxx

Bessel.o:Bessel.cxx Bessel.h Numeric.h
	${CXX} ${CXXFLAGS} -c Bessel.cxx

GeomVector.o:GeomVector.cxx GeomVector.h Numeric.h
	${CXX} ${CXXFLAGS} -c GeomVector.cxx

GeomPlane.o:GeomPlane.cxx GeomPlane.h GeomVector.h Numeric.h 
	${CXX} ${CXXFLAGS} -c GeomPlane.cxx 

IrregularGeometryInformation.o:IrregularGeometryInformation.cxx IrregularGeometryInformation.h GeomPlane.h Numeric.h \
              GeomVector.h
	${CXX} ${CXXFLAGS} -c IrregularGeometryInformation.cxx

SphericalHarmonics.o:SphericalHarmonics.cxx SphericalHarmonics.h Numeric.h
	${CXX} ${CXXFLAGS} -c SphericalHarmonics.cxx

GroupDataSet.o:GroupDataSet.cxx GroupDataSet.h GroupData.h Numeric.h
	${CXX} ${CXXFLAGS} -c GroupDataSet.cxx

Nuclide.o:Nuclide.cxx Nuclide.h GroupDataSet.h GroupData.h Numeric.h LibData.h
	${CXX} ${CXXFLAGS} -c Nuclide.cxx

DelayedNeutronData.o:DelayedNeutronData.cxx DelayedNeutronData.h Numeric.h GroupData.h Numeric.h \
	MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c DelayedNeutronData.cxx

GeneralSystem.o:GeneralSystem.cxx GeneralSystem.h Numeric.h Medium.h GeneralOption.h \
                GeneralMesh.h CartMeshInfo.h DelayedNeutronData.h DelayedNeutronData.h Nuclide.h \
		SensitivityData.h
	${CXX} ${CXXFLAGS} -c GeneralSystem.cxx

CartesianLattice.o:CartesianLattice.cxx CartesianLattice.h Numeric.h
	${CXX} ${CXXFLAGS} -c CartesianLattice.cxx

CrossSectionEdit.o:CrossSectionEdit.cxx CrossSectionEdit.h Numeric.h GeneralSystem.h Medium.h
	${CXX} ${CXXFLAGS} -c CrossSectionEdit.cxx

##

HigherModeCalculator.o:HigherModeCalculator.cxx HigherModeCalculator.h Numeric.h GeneralSystem.h \
                       DHEX_system.h PLOS_system.h

	${CXX} ${CXXFLAGS} -c HigherModeCalculator.cxx

## -- for BurnupCalculation

Burnup.o:Burnup.cxx Burnup.h Numeric.h GroupData.h GeneralSystem.h BurnupChain.h Medium.h \
	 Nuclide.h MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c Burnup.cxx

BurnupChain.o:BurnupChain.cxx BurnupChain.h Numeric.h GroupData.h GeneralSystem.h Burnup.h \
              MATIDTranslator.h BurnupChainGenerator.h
	${CXX} ${CXXFLAGS} -c BurnupChain.cxx

Burner.o:Burner.cxx Burner.h Numeric.h LibData.h Medium.h OnePointCalculator.h \
	IrregularGeometryInformation.h PJI_trajectoryset.h PJI_system.h \
	BurnupChainGenerator.h FRDesignTool.h Nuclide.h SensitivityData.h Burnup.h \
	MATIDTranslator.h SelfShieldingCalculator.h GeneralBurner.h GData.h RMforFR.h
	${CXX} ${CXXFLAGS} -c Burner.cxx

Cooler.o:Cooler.cxx Cooler.h Numeric.h ENDFIDTranslator.h BurnupChainGenerator.h Burnup.h \
	MATIDTranslator.h GroupData.h GammaTool.h GData.h
	${CXX} ${CXXFLAGS} -c Cooler.cxx

FRBurnerRZ.o:FRBurnerRZ.cxx FRBurnerRZ.h PLOSE_system.h PLOS_system.h SNRZ_system.h \
             Medium.h CartCore.h CartMeshInfo.h GeneralOption.h Burnup.h OnePointCalculator.h \
             FRDesignTool.h MATIDTranslator.h BurnupChainGenerator.h SelfShieldingCalculator.h
	${CXX} ${CXXFLAGS} -c FRBurnerRZ.cxx

## -- for library data

LibData.o:LibData.cxx LibData.h Numeric.h GroupData.h MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c LibData.cxx

## -- for Gammalibrary data

GammaLibData.o:GammaLibData.cxx GammaLibData.h Numeric.h GroupData.h MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c GammaLibData.cxx

GammaTool.o:GammaTool.cxx GammaTool.h Numeric.h GammaLibData.h Medium.h GroupData.h \
	GeneralSystem.h Burnup.h
	${CXX} ${CXXFLAGS} -c GammaTool.cxx

## -- UNC

UNC_Covariance.o:UNC_Covariance.cxx UNC_Covariance.h Numeric.h GroupData.h
	${CXX} ${CXXFLAGS} -c UNC_Covariance.cxx

UNC_CrossSection.o:UNC_CrossSection.cxx UNC_CrossSection.h UNC_Covariance.h Numeric.h GroupData.h LibData.h \
	MATIDTranslator.h 
	${CXX} ${CXXFLAGS} -c UNC_CrossSection.cxx

UNC_Parameters.o:UNC_Parameters.cxx UNC_Parameters.h UNC_Covariance.h Numeric.h GroupData.h
	${CXX} ${CXXFLAGS} -c UNC_Parameters.cxx

UNC.o:UNC.cxx UNC.h UNC_CrossSection.h UNC_Covariance.h Numeric.h MATIDTranslator.h BurnupChain.h
	${CXX} ${CXXFLAGS} -c UNC.cxx

UNC_Sensitivity.o:UNC_Sensitivity.cxx UNC_Sensitivity.h GroupData.h Numeric.h UNC_Parameters.h SensitivityData.h 
	${CXX} ${CXXFLAGS} -c UNC_Sensitivity.cxx

YieldDecayCovariance.o:YieldDecayCovariance.cxx YieldDecayCovariance.h GroupData.h Numeric.h BurnupChain.h BurnupChainGenerator.h
	${CXX} ${CXXFLAGS} -c YieldDecayCovariance.cxx

## -- for One-point calculator

OnePointCalculator.o:OnePointCalculator.cxx OnePointCalculator.h Numeric.h Medium.h LibData.h \
	Nuclide.h
	${CXX} ${CXXFLAGS} -c OnePointCalculator.cxx

SelfShieldingCalculator.o:SelfShieldingCalculator.cxx SelfShieldingCalculator.h Numeric.h Nuclide.h
	${CXX} ${CXXFLAGS} -c SelfShieldingCalculator.cxx

## -- for PLOS --

PLOS_system.o: PLOS_system.cxx PLOS_system.h GeneralSystem.h Numeric.h Medium.h SensitivityData.h
	${CXX} ${CXXFLAGS} -c PLOS_system.cxx 

## -- for PLOSE -- 

PLOSE_system.o: PLOSE_system.cxx PLOSE_system.h GeneralSystem.h Numeric.h
	${CXX} ${CXXFLAGS} -c PLOSE_system.cxx

## -- for DHEX -- 

DHEX_system.o: DHEX_system.cxx DHEX_system.h GeneralSystem.h Numeric.h SensitivityData.h
	${CXX} ${CXXFLAGS} -c DHEX_system.cxx

## -- for SNT --

SNT_system.o: SNT_system.cxx SNT_system.h Quadrature.h PLOSE_system.h Numeric.h \
              GeneralSystem.h Medium.h SensitivityData.h
	${CXX} ${CXXFLAGS} -c SNT_system.cxx

Quadrature.o: Quadrature.cxx Quadrature.h Numeric.h GroupData.h
	${CXX} ${CXXFLAGS} -c Quadrature.cxx

## -- for SNR --

SNR_system.o: SNR_system.cxx SNR_system.h SNR_quadrature.h Numeric.h GeneralSystem.h Medium.h \
              SensitivityData.h
	${CXX} ${CXXFLAGS} -c SNR_system.cxx

SNR_quadrature.o: SNR_quadrature.cxx SNR_quadrature.h Numeric.h
	${CXX} ${CXXFLAGS} -c SNR_quadrature.cxx

## -- for SNRZ --

SNRZ_system.o: SNRZ_system.cxx SNRZ_system.h SNRZ_quadrature.h PLOSE_system.h Numeric.h \
               GeneralSystem.h Medium.h SensitivityData.h
	${CXX} ${CXXFLAGS} -c SNRZ_system.cxx

SNRZ_quadrature.o: SNRZ_quadrature.cxx SNRZ_quadrature.h SNRZ_system.h Numeric.h
	${CXX} ${CXXFLAGS} -c SNRZ_quadrature.cxx

## -- for PJI --

PJI_trajectory.o: PJI_trajectory.cxx PJI_trajectory.h FunctionTable.h Numeric.h
	${CXX} ${CXXFLAGS} -c PJI_trajectory.cxx

PJI_trajectoryset.o: PJI_trajectoryset.cxx PJI_trajectoryset.h PJI_trajectory.h FunctionTable.h Numeric.h \
                     GeomVector.h GeomPlane.h IrregularGeometryInformation.h
	${CXX} ${CXXFLAGS} -c PJI_trajectoryset.cxx

PJI_slabpij.o: PJI_slabpij.cxx PJI_slabpij.h Numeric.h FunctionTable.h
	${CXX} ${CXXFLAGS} -c PJI_slabpij.cxx

PJI_system.o: PJI_system.cxx PJI_system.h  PJI_trajectoryset.h Numeric.h PJI_slabpij.h FunctionTable.h \
              GeneralSystem.h Medium.h GroupData.h Nuclide.h SensitivityData.h GData.h
	${CXX} ${CXXFLAGS} -c PJI_system.cxx

## -- for MEC --

MEC_system.o:MEC_system.cxx MEC_system.h PJI_trajectoryset.h PJI_trajectory.h \
             Numeric.h GeneralSystem.h Medium.h FunctionTable.h Quadrature.h
	${CXX} ${CXXFLAGS} -c MEC_system.cxx

## -- for PNX --

PNX_system.o:PNX_system.cxx PNX_system.h \
             Numeric.h GeneralSystem.h Medium.h FunctionTable.h
	${CXX} ${CXXFLAGS} -c PNX_system.cxx

## -- for ABEMIE --

ABEMIE_medium.o: ABEMIE_medium.cxx ABEMIE_medium.h Numeric.h
	${CXX} ${CXXFLAGS} -c ABEMIE_medium.cxx

ABEMIE_plane.o: ABEMIE_plane.cxx  ABEMIE_plane.h ABEMIE_const.h Numeric.h
	${CXX} ${CXXFLAGS} -c ABEMIE_plane.cxx

ABEMIE_region.o: ABEMIE_region.cxx ABEMIE_region.h ABEMIE_const.h ABEMIE_plane.h ABEMIE_medium.h Numeric.h
	${CXX} ${CXXFLAGS} -c ABEMIE_region.cxx

ABEMIE_core.o: ABEMIE_core.cxx ABEMIE_core.h ABEMIE_const.h ABEMIE_plane.h ABEMIE_region.h ABEMIE_medium.h Numeric.h
	${CXX} ${CXXFLAGS} -c ABEMIE_core.cxx

## -- for FRDesignTool --

FRDesignTool.o: FRDesignTool.cxx FRDesignTool.h Medium.h MATIDTranslator.h Burnup.h
	${CXX} ${CXXFLAGS} -c FRDesignTool.cxx

FRDesign.o: FRDesign.cxx FRDesign.h FRDesignTool.h Medium.h CartCore.h Burnup.h
	${CXX} ${CXXFLAGS} -c FRDesign.cxx

## -- for PromptKratio

PKR.o: PKR.cxx PKR.h Numeric.h LibData.h DelayedNeutronData.h SNR_system.h SensitivityData.h
	${CXX} ${CXXFLAGS} -c PKR.cxx

## -- for OktavianTool

OktavianTool.o: OktavianTool.cxx OktavianTool.h Numeric.h SensitivityData.h
	${CXX} ${CXXFLAGS} -c OktavianTool.cxx

## -- for SVDTool

SVDTool.o: SVDTool.cxx SVDTool.h Numeric.h GroupData.h
	${CXX} ${CXXFLAGS} -c SVDTool.cxx

## -- for ENDFIDTranslator

ENDFIDTranslator.o: ENDFIDTranslator.cxx ENDFIDTranslator.h Numeric.h
	${CXX} ${CXXFLAGS} -c ENDFIDTranslator.cxx

MATIDTranslator.o: MATIDTranslator.cxx MATIDTranslator.h Numeric.h
	${CXX} ${CXXFLAGS} -c MATIDTranslator.cxx

## --

BurnupChainGenerator.o: BurnupChainGenerator.cxx BurnupChainGenerator.h MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c BurnupChainGenerator.cxx

KineticUtility.o: KineticUtility.cxx KineticUtility.h BurnupChainGenerator.h SensitivityData.h Numeric.h \
	BurnupChain.h
	${CXX} ${CXXFLAGS} -c KineticUtility.cxx

## --

PincellCalculator.o : PincellCalculator.cxx PincellCalculator.h Numeric.h PJI_trajectoryset.h \
	Medium.h LibData.h IrregularGeometryInformation.h OnePointCalculator.h GeneralOption.h \
	SelfShieldingCalculator.h PJI_system.h SensitivityData.h MEC_system.h GData.h
	${CXX} ${CXXFLAGS} -c PincellCalculator.cxx


## --

SensitivityData.o : SensitivityData.cxx SensitivityData.h GroupData.h Numeric.h MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c SensitivityData.cxx

SUTool.o : SUTool.cxx SUTool.h GroupData.h Numeric.h SensitivityData.h LibData.h BurnupChain.h
	${CXX} ${CXXFLAGS} -c SUTool.cxx

## --

RandomSamplingSupporter.o : RandomSamplingSupporter.cxx RandomSamplingSupporter.h LibData.h \
	UNC_CrossSection.h UNC.h UNC_Covariance.h UNC_Parameters.h UNC_Sensitivity.h
	${CXX} ${CXXFLAGS} -c RandomSamplingSupporter.cxx
## --

MulticellBurner.o : MulticellBurner.cxx MulticellBurner.h \
	Numeric.h SensitivityData.h Medium.h IrregularGeometryInformation.h PJI_trajectoryset.h \
	PJI_system.h Burnup.h SelfShieldingCalculator.h OnePointCalculator.h LibData.h \
	GeneralBurner.h MEC_system.h GData.h 
	${CXX} ${CXXFLAGS} -c MulticellBurner.cxx
## --

GeneralBurner.o : GeneralBurner.cxx GeneralBurner.h \
	Numeric.h LibData.h MATIDTranslator.h OnePointCalculator.h SelfShieldingCalculator.h \
	Medium.h FRDesignTool.h PJI_trajectoryset.h Burnup.h
	${CXX} ${CXXFLAGS} -c GeneralBurner.cxx

## --

Pij_dyn.o: Pij_dyn.cxx GeomVector.h FunctionTable.h Numeric.h GeneralSystem.h BurnupChainGenerator.h
	${CXX} ${CXXFLAGS} -c Pij_dyn.cxx

## --

FMG_Tool.o: FMG_Tool.cxx FMG_Tool.h Numeric.h GroupData.h MATIDTranslator.h
	${CXX} ${CXXFLAGS} -c FMG_Tool.cxx

## -- for GeneralizedData

GData.o: GData.cxx GData.h
	${CXX} ${CXXFLAGS} -c GData.cxx

## -- for TH

HeatConductionCalculator.o: HeatConductionCalculator.cxx HeatConductionCalculator.h Numeric.h GroupData.h \
	Medium.h PLOS_system.h CartMeshInfo.h GeneralOption.h
	${CXX} ${CXXFLAGS} -c HeatConductionCalculator.cxx

jsmest.o: jsmest.f
	gfortran -c jsmest.f

jsmest_wrapper.o:jsmest_wrapper.cxx jsmest_wrapper.h 
	${CXX} ${CXXFLAGS} -c jsmest_wrapper.cxx

TwoPhaseFlowCalculator.o: TwoPhaseFlowCalculator.cxx TwoPhaseFlowCalculator.h jsmest_wrapper.h Numeric.h \
	GroupData.h
	${CXX} ${CXXFLAGS} -c TwoPhaseFlowCalculator.cxx

## -- for RMforFR

RMforFR.o: RMforFR.cxx RMforFR.h Numeric.h
	${CXX} ${CXXFLAGS} -c RMforFR.cxx

BNCTPerturbator.o: BNCTPerturbator.cxx BNCTPerturbator.h PLOS_system.h GammaTool.h Medium.h \
	SNR_system.h SNR_quadrature.cxx
	${CXX} ${CXXFLAGS} -c BNCTPerturbator.cxx 
SomeOutputMethod.o: SomeOutputMethod.cxx GeneralSystem.h GeneralMesh.h Medium.h PLOS_system.h\
	SNT_system.h CartMeshInfo.h 
	${CXX} ${CXXFLAGS} -c SomeOutputMethod.cxx 

ResonanceCalculationFan.o: ResonanceCalculationFan.cxx Medium.h Numeric.h
	${CXX} ${CXXFLAGS} -c ResonanceCalculationFan.cxx 


clean:
	${RM} ${OBJS} ${TARGET} 

distclean:
	${RM} ${OBJS} ${TARGET}

cleansolver:
	${RM} PLOS*.o SNT*.o SNR*.o SNRZ*.o 
