TOP=../..

include $(TOP)/configure/CONFIG

LIBRARY_IOC_RTEMS = mcor
LIBRARY_IOC_vxWorks = mcor

DBD += mcor.dbd

mcor_DBD += mcorSup.dbd

mcor_SRCS += mcorASub.c

mcor_LIBS += $(EPICS_BASE_IOC_LIBS)

include $(TOP)/configure/RULES
