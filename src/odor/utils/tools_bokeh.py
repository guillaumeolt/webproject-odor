from bokeh.models import ColumnDataSource, OpenURL, TapTool, Legend
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
    l_legend = list(set(l_predicted))
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
    
    
    // View of the colors for convenience
    var colors = render_cds.data['color'];
    var color_views = render_cds.data['color_view'];
    
    // Convenient to have the number of data points
    var n = colors.length;

    // Update colors
    for (var i = 0; i < n; i++) {
        if (odors[i].split(";").includes(colorby)) {
            color_views[i] = 'blue';//colors[i];
            legend_color[i] = colorby;
        }
        else if (olfactory_receptors[i].split(";").includes(colorby)) {
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

def get_bokeh_plot_odor_from_list_or(mapper, db_dict, list_or, path_svg_add=""):
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
                if or_id in chem["OlfRecept"]:
                    col = color_dict[or_id]
                    leg = or_id
            except:
                pass
        list_color.append(col)
        list_legend.append(leg)

    #BOKEH plot
    script, div = get_bokeh_plot_odor_test(mapper, list_name, list_path_svg, list_color, list_legend, list_chem_id, list_chem_odors, list_chem_or,l_predicted=list_or)
    return script, div
