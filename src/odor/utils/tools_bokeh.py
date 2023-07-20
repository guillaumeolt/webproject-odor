from bokeh.models import ColumnDataSource, OpenURL, TapTool, Legend, CustomJS, BoxAnnotation
from bokeh.plotting import figure
from bokeh.embed import components
import bokeh.plotting
import bokeh.io
#from src.Server.models import my_custom_sql_odor_get_chem, my_custom_sql, ChemicalsOdors
#from src.ServerOdors.settings import BASE_DIR
from .test import tranform_db_dict
from .tools_umap import load_umap_chem_odor


def get_bokeh_plot_odor(mapper, list_name, list_path_svg, list_color, list_legend, list_chem_id, list_chem_odors):
    """
    # create a plot
    plot = figure(plot_width=400, plot_height=400)

    # add a circle renderer with a size, color, and alpha
    plot.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)
    """

    data = dict(
        x=mapper[:, 0],
        y=mapper[:, 1],
    idChemicals = list_chem_id,
    img = list_path_svg,
    name = list_name,
    color = list_color,
    legend= list_legend,
    odors = list_chem_odors
    )

    TOOLTIPS = """
    <div>
        <div>
            Name: @name <br>
            Odors: @odors <br>
        </div>
        <div>
            <img src="@img">
        </div>

    </div>
    """


    source = ColumnDataSource(data)
    plot = figure( tooltips=TOOLTIPS)
    plot.add_tools(TapTool())
    plot.scatter('x', 'y', size=5, source=source, color='color', fill_alpha = 0.5, legend_group='legend')
    url = "../chemical/@idChemicals"
    taptool = plot.select(type=TapTool)
    taptool.callback = OpenURL(url=url)
    return components(plot)

def get_bokeh_plot_odor_test(mapper, list_name, list_path_svg, list_color, list_legend, list_chem_id, list_chem_odors, list_chem_or, l_predicted):
    """
    # create a plot
    plot = figure(plot_width=400, plot_height=400)

    # add a circle renderer with a size, color, and alpha
    plot.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)
    """

    data = dict(
        x=mapper[:, 0],
        y=mapper[:, 1],
    idChemicals = list_chem_id,
    img = list_path_svg,
    name = list_name,
    color = list_color,
    color_view=list_color,
    legend= list_legend,
    legend_color= list_legend,
    odors = list_chem_odors,
    olfactory_receptors = list_chem_or,
    )
    l_legend = list(set(list_legend))
    l_legend.append("All")
    l_legend.append("other")
    colorby_selector = bokeh.models.Select(
        title="color by", options=list(set(l_legend)), value="All", width=200,
    )

    jscode = """
    // Extract what we want to color by from selector
    var colorby = colorby_selector.value;

    // View of the legends for convenience
    var legends = render_cds.data['legend'];
    var legend_color = render_cds.data['legend_color'];
    
    
    var odors = render_cds.data['odors'];
    var olfactory_receptors = render_cds.data['olfactory_receptors'];
    var idChemicals = render_cds.data['idChemicals'];
    
    // View of the colors for convenience
    var colors = render_cds.data['color'];
    var color_views = render_cds.data['color_view'];
    
    // Convenient to have the number of data points
    var n = colors.length;

    // Update colors
    console.log(colorby);
    console.log(colorby);
    for (var i = 0; i < n; i++) {
        console.log(typeof(idChemicals[i]));
        console.log(toString(idChemicals[i]));
        if (odors[i].split(";").includes(colorby)) {
            color_views[i] = 'blue';//colors[i];
            legend_color[i] = colorby;
        }
        else if (olfactory_receptors[i].split(";").includes(colorby)) {
            color_views[i] = 'blue';//colors[i];
            legend_color[i] = colorby;
        }
        else if (idChemicals[i].toString() == colorby) {
            color_views[i] = 'blue';//colors[i];
            legend_color[i] = colorby;
        } 
        else if (colorby === 'other' && colorby === legends[i]) {
            color_views[i] = 'blue';
            legend_color[i] = "other";
        }
        else {
            color_views[i] = 'rgba( 220, 220, 220, 0.2)';
            legend_color[i] = "other";
        }
    }
    if (colorby === 'All') {
        for (var i = 0; i < n; i++) {
            color_views[i] = colors[i];
            legend_color[i] = legends[i];
        }
    }
    render_cds.change.emit();
    """
    TOOLTIPS = """
    <div>
        <div>
            Name: @name <br>
            Odors: @odors <br>
            OR: @olfactory_receptors <br>
        </div>
        <div>
            <!--@img-->
            <img src="@img">
        </div>

    </div>
    """



    render_cds = ColumnDataSource(data)
    cds = render_cds
    plot = figure( tooltips=TOOLTIPS)
    plot.add_tools(TapTool())
    plot.scatter('x', 'y', size=5, source=render_cds, color='color_view', fill_alpha = 0.5, legend_field="legend_color")
    # legend=Legend(items=[(ID[i-1],[p.renderers[0]]) for i in df.index])
    # p.add_layout(legend,'right')
    url = "../chemical/@idChemicals"
    taptool = plot.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    args = dict(
        render_cds=render_cds,
        cds=cds,
        colorby_selector=colorby_selector
    )
    colorby_selector.js_on_change("value", bokeh.models.CustomJS(code=jscode, args=args))

    layout = bokeh.layouts.row(
        plot,
        bokeh.layouts.Spacer(width=15),
        bokeh.layouts.column(
            bokeh.layouts.Spacer(height=15),
            colorby_selector,
        ),
    )

    return components(layout)

def get_bokeh_plot_odor_from_list_odors(mapper, db_dict, list_odors, path_svg_add=""):
    color_list = "orange,green,pink,blueviolet,brown,turquoise,gold,yellow,blue,red,magenta,mediumturquoise,cornflowerblue,darkblue".split(',')
    color_dict = dict(zip(list_odors, color_list[:len(list_odors)]))
    list_smi = []
    list_name = []
    list_path_svg = []
    list_color = []
    list_legend = []
    list_chem_id = []
    list_chem_odors = []
    list_chem_or = []
    for chem in db_dict:
        #with open(path_svg_add+str(chem["idChemicals"])+".svg", 'r') as file:
            #svg = file.read()
        list_name.append(chem["Name"])
        list_path_svg.append(path_svg_add+str(chem["idChemicals"])+".svg")
        list_chem_id.append(chem["idChemicals"])
        list_chem_odors.append(str(chem["smell"]))
        try:
            list_chem_or.append(str(";".join(chem["OlfRecept"])))
        except:
            list_chem_or.append("None")
        i_odor = 0
        col = 'rgba( 220, 220, 220, 0.2)'
        leg = "other"
        for odor_id in list_odors:
            try:
                if odor_id in chem["smell"].split(";"):
                    col = color_dict[odor_id]
                    leg = odor_id
            except:
                pass
        list_color.append(col)
        list_legend.append(leg)

    #BOKEH plot
    script, div = get_bokeh_plot_odor_test(mapper, list_name, list_path_svg, list_color, list_legend, list_chem_id, list_chem_odors, list_chem_or, l_predicted=list_odors)
    return script, div

def get_bokeh_plot_odor_from_list_or(mapper, db_dict, list_or, receptors_idOR_dict_bis, path_svg_add=""):
    color_list = "orange,green,pink,blueviolet,brown,turquoise,gold,yellow,blue,red,magenta,navy,skyblue,peru,cyan,darkolivegreen".split(',')
    color_dict = dict(zip(list_or, color_list[:len(list_or)]))
    list_smi = []
    list_name = []
    list_path_svg = []
    list_color = []
    list_legend = []
    list_chem_id = []
    list_chem_odors = []
    list_chem_or = []
    for chem in db_dict:
        #with open(path_svg_add+str(chem["idChemicals"])+".svg", 'r') as file:
            #svg = file.read()
        list_name.append(chem["Name"])
        list_path_svg.append(path_svg_add+str(chem["idChemicals"])+".svg")
        list_chem_id.append(chem["idChemicals"])
        list_chem_odors.append(str(chem["smell"]))
        try:
            list_chem_or.append(str(";".join(chem["OlfRecept"])))
        except:
            list_chem_or.append("None")

        i_odor = 0
        col = 'rgba( 220, 220, 220, 0.2)'
        leg = "other"
        for or_id in list_or:
            try:
                if or_id in chem["idOlfactoryReceptors"]:
                    col = color_dict[or_id]
                    leg = receptors_idOR_dict_bis[or_id]
            except:
                pass
        list_color.append(col)
        list_legend.append(leg)

    #BOKEH plot
    script, div = get_bokeh_plot_odor_test(mapper, list_name, list_path_svg, list_color, list_legend, list_chem_id, list_chem_odors, list_chem_or,l_predicted=list_or)
    return script, div



def get_bokeh_plot_odor_from_list_chem(mapper, db_dict, list_chem, path_svg_add=""):
    color_list = "orange,green,pink,blueviolet,brown,turquoise,gold,yellow,blue,red,magenta,navy,skyblue,peru,cyan,darkolivegreen".split(',')
    color_dict = dict(zip(list_chem, color_list[:len(list_chem)]))
    list_smi = []
    list_name = []
    list_path_svg = []
    list_color = []
    list_legend = []
    list_chem_id = []
    list_chem_odors = []
    list_chem_or = []
    for chem in db_dict:
        #with open(path_svg_add+str(chem["idChemicals"])+".svg", 'r') as file:
            #svg = file.read()
        list_name.append(chem["Name"])
        list_path_svg.append(path_svg_add+str(chem["idChemicals"])+".svg")
        list_chem_id.append(chem["idChemicals"])
        list_chem_odors.append(str(chem["smell"]))
        try:
            list_chem_or.append(str(";".join(chem["OlfRecept"])))
        except:
            list_chem_or.append("None")

        i_odor = 0
        col = 'rgba( 220, 220, 220, 0.2)'
        leg = "other"
        for chem_id in list_chem:
            try:
                if chem_id == str(chem["idChemicals"]):
                    col = color_dict[chem_id]
                    leg = chem_id
                    print(chem_id, col, leg , "--", chem)
            except:
                pass
        list_color.append(col)
        list_legend.append(leg)

    #BOKEH plot
    script, div = get_bokeh_plot_odor_test(mapper, list_name, list_path_svg, list_color, list_legend, list_chem_id, list_chem_odors, list_chem_or,l_predicted=list_chem)
    return script, div


def get_bokeh_plot_odor_test_bis(mapper,\
                                 list_name,\
                                 list_iupac_name,\
                                 list_pubchemcid,\
                                 list_cas,\
                                 list_smi_bis,\
                                 list_path_svg,\
                                 list_color,\
                                 list_legend,\
                                 list_chem_id,\
                                 list_chem_odors,\
                                 list_chem_or,\
                                 l_predicted):
    """
    # create a plot
    plot = figure(plot_width=400, plot_height=400)

    # add a circle renderer with a size, color, and alpha
    plot.circle([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], size=20, color="navy", alpha=0.5)
    """

    data = dict(
        x=mapper[:, 0],
        y=mapper[:, 1],
        idChemicals=list_chem_id,
        img=list_path_svg,
        name=list_name,
        iupac_name=list_iupac_name,
        pubchem_cid=list_pubchemcid,
        cas=list_cas,
        smile=list_smi_bis,
        color=list_color,
        color_view=list_color,
        legend=list_legend,
        legend_color=list_legend,
        odors=list_chem_odors,
        olfactory_receptors=list_chem_or,
    )
    l_legend = list(set(list_legend))
    l_legend.append("All")
    l_legend.append("other")
    colorby_selector = bokeh.models.Select(
        title="color by", options=list(set(l_legend)), value="All", width=200,
    )

    jscode = """
    // Extract what we want to color by from selector
    var colorby = colorby_selector.value;

    // View of the legends for convenience
    var legends = render_cds.data['legend'];
    var legend_color = render_cds.data['legend_color'];


    var odors = render_cds.data['odors'];
    var olfactory_receptors = render_cds.data['olfactory_receptors'];
    var idChemicals = render_cds.data['idChemicals'];

    // View of the colors for convenience
    var colors = render_cds.data['color'];
    var color_views = render_cds.data['color_view'];

    // Convenient to have the number of data points
    var n = colors.length;

    // Update colors
    console.log(colorby);
    console.log(colorby);
    for (var i = 0; i < n; i++) {
        console.log(typeof(idChemicals[i]));
        console.log(toString(idChemicals[i]));
        if (odors[i].split(";").includes(colorby)) {
            color_views[i] = 'blue';//colors[i];
            legend_color[i] = colorby;
        }
        else if (olfactory_receptors[i].split(";").includes(colorby)) {
            color_views[i] = 'blue';//colors[i];
            legend_color[i] = colorby;
        }
        else if (idChemicals[i].toString() == colorby) {
            color_views[i] = 'blue';//colors[i];
            legend_color[i] = colorby;
        } 
        else if (colorby === 'other' && colorby === legends[i]) {
            color_views[i] = 'blue';
            legend_color[i] = "other";
        }
        else {
            color_views[i] = 'rgba( 220, 220, 220, 0.2)';
            legend_color[i] = "other";
        }
    }
    if (colorby === 'All') {
        for (var i = 0; i < n; i++) {
            color_views[i] = colors[i];
            legend_color[i] = legends[i];
        }
    }
    render_cds.change.emit();
    """
    TOOLTIPS = """
    <style type="text/css">
    .tg  {border-collapse:collapse;border-spacing:0;border-color:black;margin-top:1em;}
    .tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:10px;
      overflow:hidden;padding:1px 3px;word-break:normal;}
    .tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:10px;
      font-weight:normal;overflow:hidden;padding:1px 3px;word-break:normal;}
    .tg .tg-l6li{border-color:#000000;font-size:10px;text-align:left;vertical-align:top}
    .tg .tg-4k6h{border-color:#000000;font-size:10px;font-weight:bold;text-align:left;vertical-align:top}
    </style>
    <table class="tg">
    <thead>
      <tr>
        <td class="tg-l6li" rowspan="7"><img src="@img" alt="molecule image" width="100" height="80"></td>
        <td class="tg-4k6h"><b>Name :</b></td>
        <td class="tg-l6li">@name</td>
      </tr>
      <tr>
        <td class="tg-4k6h"><b>IUPAC Name :</b></td>
        <td class="tg-l6li">@iupac_name</td>
      </tr>
      <tr>
        <td class="tg-4k6h"><b>PubChem CID :</b></td>
        <td class="tg-l6li">@pubchem_cid</td>
      </tr>
      <tr>
        <td class="tg-4k6h"><b>Cas N° :</b></td>
        <td class="tg-l6li">@cas</td>
      </tr>
      <tr>
        <td class="tg-4k6h"><b>SMILES :</b></td>
        <td class="tg-l6li">@smile</td>
      </tr>
      <tr>
        <td class="tg-4k6h"><b>Odors :</b></td>
        <td class="tg-l6li">@odors</td>
      </tr>
      <tr>
        <td class="tg-4k6h"><b>ORs :</b></td>
        <td class="tg-l6li">@olfactory_receptors</td>
      </tr>
    </thead>
    </table>
    
    <style type="text/css">
    .tg  {border-collapse:collapse;border-spacing:0;margin-top:1em;}
    .tg td{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:10px;
      overflow:hidden;padding:1px 3px;word-break:normal;}
    .tg th{border-color:black;border-style:solid;border-width:1px;font-family:Arial, sans-serif;font-size:10px;
      font-weight:normal;overflow:hidden;padding:1px 3px;word-break:normal;}
    .tg .tg-y3xd{border-color:#fe0000;font-size:10px;text-align:left;vertical-align:middle}
    .tg .tg-magg{border-color:#fe0000;font-size:10px;font-weight:bold;text-align:left;vertical-align:middle}
    </style>
    """

    render_cds = ColumnDataSource(data)
    cds = render_cds
    plot = figure(title='Pan and Zoom Here', tooltips=TOOLTIPS)
    plot.add_tools(TapTool())
    plot.scatter('x', 'y', size=5, source=render_cds, color='color_view', fill_alpha=0.5, legend_field="legend_color")
    # legend=Legend(items=[(ID[i-1],[p.renderers[0]]) for i in df.index])
    # p.add_layout(legend,'right')
    url = "../chemical/@idChemicals"
    taptool = plot.select(type=TapTool)
    taptool.callback = OpenURL(url=url)

    # box 2
    box = BoxAnnotation(left=0, right=0, bottom=0, top=0,
                        fill_alpha=0.1, line_color='black', fill_color='black')

    jscode2 = """
        box[%r] = cb_obj.start
        box[%r] = cb_obj.end
    """
    xcb = CustomJS(args=dict(box=box), code=jscode2 % ('left', 'right'))
    ycb = CustomJS(args=dict(box=box), code=jscode2 % ('bottom', 'top'))
    plot.x_range.js_on_change('start', xcb)
    plot.x_range.js_on_change('end', xcb)
    plot.y_range.js_on_change('start', ycb)
    plot.y_range.js_on_change('end', ycb)

    p2 = figure(title='See Zoom Window Here', tools='')
    p2.scatter('x', 'y', size=5, source=render_cds, color='color_view', fill_alpha=0.5, legend_field="legend_color")
    p2.add_layout(box)

    args = dict(
        render_cds=render_cds,
        cds=cds,
        colorby_selector=colorby_selector
    )
    colorby_selector.js_on_change("value", bokeh.models.CustomJS(code=jscode, args=args))

    layout = bokeh.layouts.row(
        plot,
        bokeh.layouts.Spacer(width=15),
        bokeh.layouts.column(
            bokeh.layouts.Spacer(height=15),
            colorby_selector,
        ),
        p2
    )

    return components(layout)

def get_bokeh_plot_odor_from_list_or_bis(mapper, db_dict, list_or, receptors_idOR_dict_bis, path_svg_add=""):
    color_list = "orange,green,pink,blueviolet,brown,turquoise,gold,yellow,blue,red,magenta,navy,skyblue,peru,cyan,darkolivegreen".split(',')
    color_dict = dict(zip(list_or, color_list[:len(list_or)]))
    list_smi = []
    list_name = []
    list_smi_bis = []
    list_name = []
    list_iupac_name = []
    list_pubchemcid = []
    list_cas = []
    list_path_svg = []
    list_color = []
    list_legend = []
    list_chem_id = []
    list_chem_odors = []
    list_chem_or = []
    for chem in db_dict:
        #with open(path_svg_add+str(chem["idChemicals"])+".svg", 'r') as file:
            #svg = file.read()
        list_name.append(chem["Name"])
        list_iupac_name.append(chem["IUPAC_name"])
        list_pubchemcid.append(chem["Pubchem_CID"])
        list_cas.append(chem["CAS"])
        list_smi_bis.append(chem["SMILE"])
        list_path_svg.append(path_svg_add+str(chem["idChemicals"])+".svg")
        list_chem_id.append(chem["idChemicals"])
        list_chem_odors.append(str(chem["smell"]))
        try:
            list_chem_or.append(str(";".join(chem["OlfRecept"])))
        except:
            list_chem_or.append("None")

        i_odor = 0
        col = 'rgba( 220, 220, 220, 0.2)'
        leg = "other"
        for or_id in list_or:
            try:
                if or_id in chem["idOlfactoryReceptors"]:
                    col = color_dict[or_id]
                    leg = receptors_idOR_dict_bis[or_id]
            except:
                pass
        list_color.append(col)
        list_legend.append(leg)

    #BOKEH plot
    script, div = get_bokeh_plot_odor_test_bis(mapper,\
                                               list_name,\
                                               list_iupac_name,\
                                               list_pubchemcid,\
                                               list_cas,\
                                               list_smi_bis,\
                                               list_path_svg,\
                                               list_color,\
                                               list_legend,\
                                               list_chem_id,\
                                               list_chem_odors,\
                                               list_chem_or,l_predicted=list_or)
    return script, div

def get_bokeh_plot_odor_from_list_chem_bis(mapper, db_dict, list_chem, path_svg_add=""):
    color_list = "orange,green,pink,blueviolet,brown,turquoise,gold,yellow,blue,red,magenta,navy,skyblue,peru,cyan,darkolivegreen".split(',')
    color_dict = dict(zip(list_chem, color_list[:len(list_chem)]))
    list_smi = []
    list_smi_bis = []
    list_name = []
    list_iupac_name = []
    list_pubchemcid = []
    list_cas = []
    list_path_svg = []
    list_color = []
    list_legend = []
    list_chem_id = []
    list_chem_odors = []
    list_chem_or = []
    for chem in db_dict:
        #with open(path_svg_add+str(chem["idChemicals"])+".svg", 'r') as file:
            #svg = file.read()
        list_name.append(chem["Name"])
        list_iupac_name.append(chem["IUPAC_name"])
        list_pubchemcid.append(chem["Pubchem_CID"])
        list_cas.append(chem["CAS"])
        list_smi_bis.append(chem["SMILE"])

        list_path_svg.append(path_svg_add+str(chem["idChemicals"])+".svg")
        list_chem_id.append(chem["idChemicals"])
        list_chem_odors.append(str(chem["smell"]))
        try:
            list_chem_or.append(str(";".join(chem["OlfRecept"])))
        except:
            list_chem_or.append("None")

        i_odor = 0
        col = 'rgba( 220, 220, 220, 0.2)'
        leg = "other"
        for chem_id in list_chem:
            try:
                if chem_id == str(chem["idChemicals"]):
                    col = color_dict[chem_id]
                    leg = chem_id
                    #print(chem_id, col, leg , "--", chem)
            except:
                pass
        list_color.append(col)
        list_legend.append(leg)

    #BOKEH plot
    script, div = get_bokeh_plot_odor_test_bis(mapper,\
                                               list_name,\
                                               list_iupac_name,\
                                               list_pubchemcid,\
                                               list_cas,\
                                               list_smi_bis,\
                                               list_path_svg,\
                                               list_color,\
                                               list_legend,\
                                               list_chem_id,\
                                               list_chem_odors,\
                                               list_chem_or,l_predicted=list_chem)
    return script, div

def get_bokeh_plot_odor_from_list_odors_bis(mapper, db_dict, list_odors, path_svg_add=""):
    color_list = "orange,green,pink,blueviolet,brown,turquoise,gold,yellow,blue,red,magenta,mediumturquoise,cornflowerblue,darkblue".split(',')
    color_dict = dict(zip(list_odors, color_list[:len(list_odors)]))
    list_smi = []
    list_name = []
    list_smi_bis = []
    list_name = []
    list_iupac_name = []
    list_pubchemcid = []
    list_cas = []
    list_path_svg = []
    list_color = []
    list_legend = []
    list_chem_id = []
    list_chem_odors = []
    list_chem_or = []
    for chem in db_dict:
        #with open(path_svg_add+str(chem["idChemicals"])+".svg", 'r') as file:
            #svg = file.read()
        list_name.append(chem["Name"])
        list_iupac_name.append(chem["IUPAC_name"])
        list_pubchemcid.append(chem["Pubchem_CID"])
        list_cas.append(chem["CAS"])
        list_smi_bis.append(chem["SMILE"])
        list_path_svg.append(path_svg_add+str(chem["idChemicals"])+".svg")
        list_chem_id.append(chem["idChemicals"])
        list_chem_odors.append(str(chem["smell"]))
        try:
            list_chem_or.append(str(";".join(chem["OlfRecept"])))
        except:
            list_chem_or.append("None")
        i_odor = 0
        col = 'rgba( 220, 220, 220, 0.2)'
        leg = "other"
        for odor_id in list_odors:
            try:
                if odor_id in chem["smell"].split(";"):
                    col = color_dict[odor_id]
                    leg = odor_id
            except:
                pass
        list_color.append(col)
        list_legend.append(leg)

    #BOKEH plot
    script, div = get_bokeh_plot_odor_test_bis(mapper,
                                               list_name,
                                               list_iupac_name, \
                                               list_pubchemcid, \
                                               list_cas, \
                                               list_smi_bis, \
                                               list_path_svg,
                                               list_color,
                                               list_legend,
                                               list_chem_id,
                                               list_chem_odors,
                                               list_chem_or, l_predicted=list_odors)
    return script, div