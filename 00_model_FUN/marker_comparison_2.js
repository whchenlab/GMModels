/**
 * =======================================================================================
 *   last modified: Jan 15, 2021; 
 *  + 由调用函数控制所有全局变量
 * =======================================================================================
 */
// 获取输入数据 ...
var alldata = JSON.parse(data);
var data = alldata.data; // 真正的数据

/**
 * =======================================================================================
 *   全局变量 ... 可根据要求进行调整，特别是 nproj_cutoff 和 lda_cutoff
 * =======================================================================================
 */
var lPatchWidth = 200;
var itemSize = 12,
    cellMargin = 0.5,
    cellSize = itemSize - cellMargin * 2,
    margin = {
        top: 20,
        right: 20,
        bottom: 320,
        left: 130
    },
    axisTickLabelFontSize = alldata.axisTickLabelFontSize, // xy轴字体大小
    lda_cutoff = alldata.lda_cutoff, // lda cutoff? 默认为 2 
    nproj_cutoff = alldata.nproj_cutoff; // 一个marker在多少个项目当中出现？
    
var colorHold = ["#2E6E12", "#57C84D", "#83D475", "#C5E8B7", "#F6BDC0", "#F1959B", "#EA4C46", "#781426"];
var colorLText = ["< -4", "-3", "-2", "-1", "1", "2", "3", "> 4"];

/**
 * =======================================================================================
 *   子函数 ...
 * =======================================================================================
 */
function bandClassifier(val) {
    var newval = val < 0 ? Math.floor(val) + 2 + 3 : Math.floor(val) + 3;
    var newval2 = newval < 0 ? 0 : newval > 7 ? 7 : newval;
    return newval2;
}



/**
 * =======================================================================================
 *    获取XY 轴；
 *    X：marker taxa;
 *    Y: projects
 *  Nov 21, 2020: add filter according to nrproj --
 * =======================================================================================
 */
var markertaxa = d3.set(data.filter(function(d) {
        return (d.nrproj >= nproj_cutoff) && (Math.abs(d.LDA) >= lda_cutoff)
    }).map(function(item) {
        return item.scientific_name;
    })).values();
var projects = d3.set(data.map(function(item) {
        return item.project_id;
    })).values();

/**
 * =======================================================================================
 *   决定画图的长宽
 * =======================================================================================
 */
var width = itemSize * markertaxa.length + 500,
    height = itemSize * projects.length;

svg
    .attr("width", width + margin.left + margin.right)
    .attr("height", height + margin.top + margin.bottom)
    .attr("font-family", "sans-serif");

/**
 * #################################################################
 *      SVG contents ...
 * #################################################################
 */
var canvas = svg.append("g")
    .attr("transform", "translate(" + margin.left + "," + margin.top + ")");

/**
 * --------------------------
 * y axis, manual plot
 * --------------------------
 */
var yaxisG = canvas.append("g").attr("id", "yaxis");
var yaxisLabel = yaxisG.selectAll("yaxislabels").data(projects).enter().append("g");

yaxisLabel.append("text")
    .attr('font-size', axisTickLabelFontSize)
    .attr("x", -5)
    .attr("y", function(d, i) {
        return itemSize * (i + 0.5) + axisTickLabelFontSize * 1 / 3;
    })
    .attr("text-anchor", "end")
    .text(function(d) {
        return d;
    });
    
yaxisLabel.style("cursor", "pointer");
yaxisLabel.append("line")
    .attr("x1", -3)
    .attr("y1", function(d, i) {
        return itemSize * (i + 0.5);
    })
    .attr("x2", 0)
    .attr("y2", function(d, i) {
        return itemSize * (i + 0.5);
    })
    .style("stroke", "black");

/**
 * --------------------------
 * x axis, manual plot
 * --------------------------
 */

var xaxisG = canvas.append("g").attr("id", "xaxis");
var xaxisLabel = xaxisG.selectAll("xaxislabels").data(markertaxa).enter().append("g");

xaxisLabel.append("text")
    .style("text-anchor", "end")
    .attr('font-size', axisTickLabelFontSize)
    .attr("transform", function(d, i) {
        var x = itemSize * (i + 0.5) + axisTickLabelFontSize * 1 / 3;
        var y = projects.length * itemSize + 8;
        return "translate(" + x + "," + y + ") rotate(-65)";
    })
    .text(function(d) {
        return d;
    });

xaxisLabel.append("line")
    .attr("x1", function(d, i) {
        return itemSize * (i + 0.5);
    })
    .attr("y1", projects.length * itemSize + 3)
    .attr("x2", function(d, i) {
        return itemSize * (i + 0.5);
    })
    .attr("y2", projects.length * itemSize)
    .style("stroke", "black");

/**
 * --------------------------
 * plot TILEs
 *  Nov 21, 2020: add filter according nr of projects in which the taxa are marker
 * --------------------------
 */
var cells = canvas.selectAll('rect')
    .data(data.filter(function(d) {
        return d.nrproj >= nproj_cutoff && Math.abs(d.LDA) >= lda_cutoff
    }))
    .enter().append('g').append('rect')
    .attr('class', 'cell')
    .attr('width', cellSize)
    .attr('height', cellSize)
    .attr('y', function(d) {
        return projects.indexOf(d.project_id) * itemSize + cellMargin;
    })
    .attr('x', function(d) {
        return markertaxa.indexOf(d.scientific_name) * itemSize + cellMargin;
    })
    .attr('fill', function(d) {
        return colorHold[bandClassifier(d.LDA)];
    })
    .attr("rx", 3)
    .attr("ry", 3)
    .attr("title", function(d) {
        return d.LDA
    });

/**
 * --------------------------
 * legend
 * --------------------------
 */
var legends = canvas.append("g").attr("class", "legends")
    .attr("transform", "translate(" + 0 + "," + (height + margin.bottom - 35 - 20) + ")");

var legend_tile_width = lPatchWidth / colorLText.length;

// title ...
legends.append("text")
    .attr("x", lPatchWidth / 2)
    .attr('font-weight', 'bold')
    .attr('font-size', axisTickLabelFontSize)
    .style("text-anchor", "middle")
    .text("LDA score")

//Legend Rectangels
legends.append("g").attr("class", "LegRect")
    .attr("transform", "translate(0," + 10 + ")")
    .selectAll("rect").data(colorHold).enter()
    .append("rect").attr("width", lPatchWidth / colorHold.length + "px").attr("height", "10px").attr("fill", function(d) {
        return d
    })
    .attr("x", function(d, i) {
        return i * (lPatchWidth / colorHold.length)
    });

// legend text
legends.append("g").attr("class", "LegText")
    .attr("transform", "translate(0," + 10 + ")")
    .selectAll("text").data(colorLText).enter()
    .append("text")
    .attr("x", function(d, i) {
        return (i + 0.5) * legend_tile_width
    })
    .attr("y", axisTickLabelFontSize + 10)
    .attr('font-weight', 'normal')
    .style("text-anchor", "middle")
    .text(function(d) {
        return d;
    });
