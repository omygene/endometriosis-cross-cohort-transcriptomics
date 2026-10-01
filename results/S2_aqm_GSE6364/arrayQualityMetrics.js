// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = ["false","true","false","false","false","false","false","false","false","true","false","false","false","false","false","true","false","false","false","false","false","false","false","true","false","true","false","false","false","false","false","false","false","false","false","false","false"];
var arrayMetadata    = [{"array":"1","sampleNames":"GSM150190.CEL.gz","index":"1"},{"array":"2","sampleNames":"GSM150191.CEL.gz","index":"2"},{"array":"3","sampleNames":"GSM150192.CEL.gz","index":"3"},{"array":"4","sampleNames":"GSM150193.CEL.gz","index":"4"},{"array":"5","sampleNames":"GSM150194.CEL.gz","index":"5"},{"array":"6","sampleNames":"GSM150195.CEL.gz","index":"6"},{"array":"7","sampleNames":"GSM150196.CEL.gz","index":"7"},{"array":"8","sampleNames":"GSM150197.CEL.gz","index":"8"},{"array":"9","sampleNames":"GSM150198.CEL.gz","index":"9"},{"array":"10","sampleNames":"GSM150199.CEL.gz","index":"10"},{"array":"11","sampleNames":"GSM150201.CEL.gz","index":"11"},{"array":"12","sampleNames":"GSM150202.CEL.gz","index":"12"},{"array":"13","sampleNames":"GSM150203.CEL.gz","index":"13"},{"array":"14","sampleNames":"GSM150204.CEL.gz","index":"14"},{"array":"15","sampleNames":"GSM150205.CEL.gz","index":"15"},{"array":"16","sampleNames":"GSM150206.CEL.gz","index":"16"},{"array":"17","sampleNames":"GSM150207.CEL.gz","index":"17"},{"array":"18","sampleNames":"GSM150208.CEL.gz","index":"18"},{"array":"19","sampleNames":"GSM150209.CEL.gz","index":"19"},{"array":"20","sampleNames":"GSM150210.CEL.gz","index":"20"},{"array":"21","sampleNames":"GSM150211.CEL.gz","index":"21"},{"array":"22","sampleNames":"GSM150212.CEL.gz","index":"22"},{"array":"23","sampleNames":"GSM150213.CEL.gz","index":"23"},{"array":"24","sampleNames":"GSM150214.CEL.gz","index":"24"},{"array":"25","sampleNames":"GSM150215.CEL.gz","index":"25"},{"array":"26","sampleNames":"GSM150216.CEL.gz","index":"26"},{"array":"27","sampleNames":"GSM150217.CEL.gz","index":"27"},{"array":"28","sampleNames":"GSM150218.CEL.gz","index":"28"},{"array":"29","sampleNames":"GSM150219.CEL.gz","index":"29"},{"array":"30","sampleNames":"GSM150220.CEL.gz","index":"30"},{"array":"31","sampleNames":"GSM150221.CEL.gz","index":"31"},{"array":"32","sampleNames":"GSM150222.CEL.gz","index":"32"},{"array":"33","sampleNames":"GSM150223.CEL.gz","index":"33"},{"array":"34","sampleNames":"GSM150224.CEL.gz","index":"34"},{"array":"35","sampleNames":"GSM150225.CEL.gz","index":"35"},{"array":"36","sampleNames":"GSM150226.CEL.gz","index":"36"},{"array":"37","sampleNames":"GSM150227.CEL.gz","index":"37"}];
var svgObjectNames   = ["pca","dens"];

var cssText = ["stroke-width:1; stroke-opacity:0.4",
               "stroke-width:3; stroke-opacity:1" ];

// Global variables - these are set up below by 'reportinit'
var tables;             // array of all the associated ('tooltips') tables on the page
var checkboxes;         // the checkboxes
var ssrules;


function reportinit() 
{
 
    var a, i, status;

    /*--------find checkboxes and set them to start values------*/
    checkboxes = document.getElementsByName("ReportObjectCheckBoxes");
    if(checkboxes.length != highlightInitial.length)
	throw new Error("checkboxes.length=" + checkboxes.length + "  !=  "
                        + " highlightInitial.length="+ highlightInitial.length);
    
    /*--------find associated tables and cache their locations------*/
    tables = new Array(svgObjectNames.length);
    for(i=0; i<tables.length; i++) 
    {
        tables[i] = safeGetElementById("Tab:"+svgObjectNames[i]);
    }

    /*------- style sheet rules ---------*/
    var ss = document.styleSheets[0];
    ssrules = ss.cssRules ? ss.cssRules : ss.rules; 

    /*------- checkboxes[a] is (expected to be) of class HTMLInputElement ---*/
    for(a=0; a<checkboxes.length; a++)
    {
	checkboxes[a].checked = highlightInitial[a];
        status = checkboxes[a].checked; 
        setReportObj(a+1, status, false);
    }

}


function safeGetElementById(id)
{
    res = document.getElementById(id);
    if(res == null)
        throw new Error("Id '"+ id + "' not found.");
    return(res)
}

/*------------------------------------------------------------
   Highlighting of Report Objects 
 ---------------------------------------------------------------*/
function setReportObj(reportObjId, status, doTable)
{
    var i, j, plotObjIds, selector;

    if(doTable) {
	for(i=0; i<svgObjectNames.length; i++) {
	    showTipTable(i, reportObjId);
	} 
    }

    /* This works in Chrome 10, ssrules will be null; we use getElementsByClassName and loop over them */
    if(ssrules == null) {
	elements = document.getElementsByClassName("aqm" + reportObjId); 
	for(i=0; i<elements.length; i++) {
	    elements[i].style.cssText = cssText[0+status];
	}
    } else {
    /* This works in Firefox 4 */
    for(i=0; i<ssrules.length; i++) {
        if (ssrules[i].selectorText == (".aqm" + reportObjId)) {
		ssrules[i].style.cssText = cssText[0+status];
		break;
	    }
	}
    }

}

/*------------------------------------------------------------
   Display of the Metadata Table
  ------------------------------------------------------------*/
function showTipTable(tableIndex, reportObjId)
{
    var rows = tables[tableIndex].rows;
    var a = reportObjId - 1;

    if(rows.length != arrayMetadata[a].length)
	throw new Error("rows.length=" + rows.length+"  !=  arrayMetadata[array].length=" + arrayMetadata[a].length);

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = arrayMetadata[a][i];
}

function hideTipTable(tableIndex)
{
    var rows = tables[tableIndex].rows;

    for(i=0; i<rows.length; i++) 
 	rows[i].cells[1].innerHTML = "";
}


/*------------------------------------------------------------
  From module 'name' (e.g. 'density'), find numeric index in the 
  'svgObjectNames' array.
  ------------------------------------------------------------*/
function getIndexFromName(name) 
{
    var i;
    for(i=0; i<svgObjectNames.length; i++)
        if(svgObjectNames[i] == name)
	    return i;

    throw new Error("Did not find '" + name + "'.");
}


/*------------------------------------------------------------
  SVG plot object callbacks
  ------------------------------------------------------------*/
function plotObjRespond(what, reportObjId, name)
{

    var a, i, status;

    switch(what) {
    case "show":
	i = getIndexFromName(name);
	showTipTable(i, reportObjId);
	break;
    case "hide":
	i = getIndexFromName(name);
	hideTipTable(i);
	break;
    case "click":
        a = reportObjId - 1;
	status = !checkboxes[a].checked;
	checkboxes[a].checked = status;
	setReportObj(reportObjId, status, true);
	break;
    default:
	throw new Error("Invalid 'what': "+what)
    }
}

/*------------------------------------------------------------
  checkboxes 'onchange' event
------------------------------------------------------------*/
function checkboxEvent(reportObjId)
{
    var a = reportObjId - 1;
    var status = checkboxes[a].checked;
    setReportObj(reportObjId, status, true);
}


/*------------------------------------------------------------
  toggle visibility
------------------------------------------------------------*/
function toggle(id){
  var head = safeGetElementById(id + "-h");
  var body = safeGetElementById(id + "-b");
  var hdtxt = head.innerHTML;
  var dsp;
  switch(body.style.display){
    case 'none':
      dsp = 'block';
      hdtxt = '-' + hdtxt.substr(1);
      break;
    case 'block':
      dsp = 'none';
      hdtxt = '+' + hdtxt.substr(1);
      break;
  }  
  body.style.display = dsp;
  head.innerHTML = hdtxt;
}
