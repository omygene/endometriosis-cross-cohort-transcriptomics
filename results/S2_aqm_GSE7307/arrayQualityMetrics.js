// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = ["false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","true","false","false","false","false","false","false","false","false","false","true","false","true","false"];
var arrayMetadata    = [{"array":"1","sampleNames":"GSM175786.CEL.gz","index":"1","group":"endometriosis"},{"array":"2","sampleNames":"GSM175787.CEL.gz","index":"2","group":"endometriosis"},{"array":"3","sampleNames":"GSM175788.CEL.gz","index":"3","group":"endometriosis"},{"array":"4","sampleNames":"GSM176039.CEL.gz","index":"4","group":"normal_endometrium"},{"array":"5","sampleNames":"GSM176040.CEL.gz","index":"5","group":"normal_endometrium"},{"array":"6","sampleNames":"GSM176041.CEL.gz","index":"6","group":"normal_endometrium"},{"array":"7","sampleNames":"GSM176043.CEL.gz","index":"7","group":"normal_endometrium"},{"array":"8","sampleNames":"GSM176082.CEL.gz","index":"8","group":"endometriosis"},{"array":"9","sampleNames":"GSM176083.CEL.gz","index":"9","group":"endometriosis"},{"array":"10","sampleNames":"GSM176084.CEL.gz","index":"10","group":"endometriosis"},{"array":"11","sampleNames":"GSM176085.CEL.gz","index":"11","group":"endometriosis"},{"array":"12","sampleNames":"GSM176086.CEL.gz","index":"12","group":"endometriosis"},{"array":"13","sampleNames":"GSM176087.CEL.gz","index":"13","group":"endometriosis"},{"array":"14","sampleNames":"GSM176088.CEL.gz","index":"14","group":"endometriosis"},{"array":"15","sampleNames":"GSM176089.CEL.gz","index":"15","group":"endometriosis"},{"array":"16","sampleNames":"GSM176090.CEL.gz","index":"16","group":"endometriosis"},{"array":"17","sampleNames":"GSM176091.CEL.gz","index":"17","group":"endometriosis"},{"array":"18","sampleNames":"GSM176092.CEL.gz","index":"18","group":"endometriosis"},{"array":"19","sampleNames":"GSM176093.CEL.gz","index":"19","group":"normal_endometrium"},{"array":"20","sampleNames":"GSM176094.CEL.gz","index":"20","group":"normal_endometrium"},{"array":"21","sampleNames":"GSM176095.CEL.gz","index":"21","group":"normal_endometrium"},{"array":"22","sampleNames":"GSM176096.CEL.gz","index":"22","group":"normal_endometrium"},{"array":"23","sampleNames":"GSM176097.CEL.gz","index":"23","group":"normal_endometrium"},{"array":"24","sampleNames":"GSM176098.CEL.gz","index":"24","group":"normal_endometrium"},{"array":"25","sampleNames":"GSM176099.CEL.gz","index":"25","group":"normal_endometrium"},{"array":"26","sampleNames":"GSM176100.CEL.gz","index":"26","group":"normal_endometrium"},{"array":"27","sampleNames":"GSM176101.CEL.gz","index":"27","group":"normal_endometrium"},{"array":"28","sampleNames":"GSM176127.CEL.gz","index":"28","group":"normal_endometrium"},{"array":"29","sampleNames":"GSM176132.CEL.gz","index":"29","group":"normal_endometrium"},{"array":"30","sampleNames":"GSM176137.CEL.gz","index":"30","group":"normal_endometrium"},{"array":"31","sampleNames":"GSM176141.CEL.gz","index":"31","group":"normal_endometrium"},{"array":"32","sampleNames":"GSM176142.CEL.gz","index":"32","group":"normal_endometrium"},{"array":"33","sampleNames":"GSM176143.CEL.gz","index":"33","group":"normal_endometrium"},{"array":"34","sampleNames":"GSM176144.CEL.gz","index":"34","group":"normal_endometrium"},{"array":"35","sampleNames":"GSM176145.CEL.gz","index":"35","group":"normal_endometrium"},{"array":"36","sampleNames":"GSM176146.CEL.gz","index":"36","group":"normal_endometrium"},{"array":"37","sampleNames":"GSM176234.CEL.gz","index":"37","group":"endometriosis"},{"array":"38","sampleNames":"GSM176236.CEL.gz","index":"38","group":"endometriosis"},{"array":"39","sampleNames":"GSM176238.CEL.gz","index":"39","group":"endometriosis"},{"array":"40","sampleNames":"GSM176240.CEL.gz","index":"40","group":"endometriosis"},{"array":"41","sampleNames":"GSM176319.CEL.gz","index":"41","group":"normal_endometrium"}];
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
