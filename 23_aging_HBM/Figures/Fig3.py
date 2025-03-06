import pandas as pd

from pyecharts import options as opts
from pyecharts.charts import Sankey

color1 = '#BAE0CA'
color2 = '#75B29A'
color3 = '#47BBC4'
color4 = '#E6DEC7'


sankey=(
        Sankey()
        .add(series_name='毕业生流向' ###给个桑基宝宝取个名字
             ,nodes=[{'name':'Frontal lobe','itemStyle':{"color": '#C8C8C8'}}
                     ,{'name':'Temporal lobe', 'itemStyle':{"color": '#C8C8C8'}}
                     ,{'name':'Parietal lobe', 'itemStyle':{"color": '#C8C8C8'}}
                     ,{'name':'Insular lobe', 'itemStyle':{"color": '#C8C8C8'}}
                     ,{'name':'Limbic lobe', 'itemStyle':{"color": '#C8C8C8'}}
                     ,{'name':'Occipital lobe', 'itemStyle':{"color": '#C8C8C8'}}
                     ,{'name':'Subcortical nuclei', 'itemStyle':{"color": '#C8C8C8'}}
                     ,{'name':'1-1', 'itemStyle':{"color": color1}}
                     ,{'name':'1-2', 'itemStyle':{"color": color2}}
                     ,{'name':'1-3', 'itemStyle':{"color": color3}}
                     ,{'name':'1-4', 'itemStyle':{"color": color4}}
                     ,{'name':'2-1', 'itemStyle':{"color": color1}}
                     ,{'name':'2-2', 'itemStyle':{"color": color2}}
                     ,{'name':'2-3', 'itemStyle':{"color": color3}}
                     ,{'name':'2-4', 'itemStyle':{"color": color4}}
                     ,{'name':'3-1', 'itemStyle':{"color": color1}} 
                     ,{'name':'3-2', 'itemStyle':{"color": color2}} 
                     ,{'name':'3-3', 'itemStyle':{"color": color3}} 
                     ,{'name':'3-4', 'itemStyle':{"color": color4}} 
                     
                      ]   ##配置有多少个节点
             ,links=[
                    {'source':'Frontal lobe','target':'1-1','value':381}
                    ,{'source':'Frontal lobe','target':'1-2','value':13855}
                    ,{'source':'Frontal lobe','target':'1-3','value':3240}
                    ,{'source':'Temporal lobe','target':'1-1','value':1681}
                    ,{'source':'Temporal lobe','target':'1-2','value':8811}
                    ,{'source':'Temporal lobe','target':'1-3','value':3781}
                    ,{'source':'Temporal lobe','target':'1-4','value':119}
                    ,{'source':'Parietal lobe','target':'1-2','value':5067}
                    ,{'source':'Parietal lobe','target':'1-3','value':4688}
                    ,{'source':'Parietal lobe','target':'1-4','value':11}
                    ,{'source':'Insular lobe','target':'1-1','value':101}
                    ,{'source':'Insular lobe','target':'1-2','value':2622}
                    ,{'source':'Insular lobe','target':'1-3','value':361}
                    ,{'source':'Limbic lobe','target':'1-2','value':2920}
                    ,{'source':'Limbic lobe','target':'1-3','value':665}
                    ,{'source':'Limbic lobe','target':'1-4','value':13}
                    ,{'source':'Occipital lobe','target':'1-2','value':1127}
                    ,{'source':'Occipital lobe','target':'1-3','value':4527}
                    ,{'source':'Subcortical nuclei','target':'1-1','value':1882}
                    ,{'source':'Subcortical nuclei','target':'1-2','value':6029}
                    ,{'source':'Subcortical nuclei','target':'1-3','value':1341}
                    ,{'source':'1-1','target':'2-1','value':3145}
                    ,{'source':'1-1','target':'2-2','value':900}
                    ,{'source':'1-2','target':'2-1','value':575}
                    ,{'source':'1-2','target':'2-2','value':36637}
                    ,{'source':'1-2','target':'2-3','value':3219}
                    ,{'source':'1-3','target':'2-2','value':2709}
                    ,{'source':'1-3','target':'2-3','value':15892}
                    ,{'source':'1-4','target':'2-3','value':65}
                    ,{'source':'1-4','target':'2-4','value':78}
                    ,{'source':'2-1','target':'3-1','value':3077}
                    ,{'source':'2-1','target':'3-2','value':643}
                    ,{'source':'2-2','target':'3-1','value':828}
                    ,{'source':'2-2','target':'3-2','value':38545}
                    ,{'source':'2-2','target':'3-3','value':873}
                    ,{'source':'2-3','target':'3-2','value':5080}
                    ,{'source':'2-3','target':'3-3','value':14081}
                    ,{'source':'2-3','target':'3-4','value':15}
                    ,{'source':'2-4','target':'3-4','value':80}


                     ]   ###配置节点之间的信息流关系
             ,linestyle_opt=opts.LineStyleOpts(opacity=0.5 ###透明度设置
                                               , curve=0.5  ###信息流的曲线弯曲度设置
                                               ,color="source" ##颜色设置，source表示使用节点的颜色
                                               ) ##线条格式 ,设置所有线条的格式
             ,label_opts=opts.LabelOpts(font_size=16
                                        ,position='left'
                                        ) ##标签配置，具体参数详见opts.LabelOpts()

           ,levels= [opts.SankeyLevelsOpts(  depth=0, ##第一层的配置
                                            
                                            linestyle_opts=opts.LineStyleOpts(color="target", opacity=0.5, curve=0.5))##信息流的配置
                    ,opts.SankeyLevelsOpts(  depth=1,##第二层的配置
                                            
                                            linestyle_opts=opts.LineStyleOpts(color="target", opacity=0.5, curve=0.5))##信息的配置
                    ,opts.SankeyLevelsOpts(  depth=2,##第三层的配置
                                            
                                            linestyle_opts=opts.LineStyleOpts(color="target", opacity=0.5, curve=0.5))##信息的配置
                    ,opts.SankeyLevelsOpts(  depth=3,##第四层的配置
                                           
                                            linestyle_opts=opts.LineStyleOpts(color="target", opacity=0.5, curve=0.5))##信息的配置
      
                    
                   ]# 桑基图每一层的设置。可以逐层设置
                )
        )
sankey.render('E:\GS Coherence\FIG_revision\Fig2\Sankey.html')