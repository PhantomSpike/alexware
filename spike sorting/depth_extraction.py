from phycontrib.template import TemplateController
from phy.utils._misc import _read_python
tc = TemplateController(**_read_python('/media/phant0msp1ke/4TB SSD/DATA/Ferret_exp_17_10_2017/Penetration2/P02-fine_quning.6/CRA/P02-fine_quning.6_cra_allch_3000Hz_768clusters/params.py'))
#for ii in range(0, 100000):
#    try:
#        print(tc.get_best_channel(ii))  # print the best channel for cluster #10
#    except:
#        break
#tc.model.get_template(10)  # gives an object with the best channel, the template waveform...