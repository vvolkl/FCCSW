from resmaster import *
pgun.energyMin = 100000
pgun.energyMax = 100000
import uuid
out.filename = "output_fullCalo_SimAndDigi_e"+ str(pgun.energyMin / 1000) +"GeV_1000events_"+uuid.uuid4().hex+".root"
