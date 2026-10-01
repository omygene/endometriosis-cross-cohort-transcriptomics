// (C) Wolfgang Huber 2010-2011

// Script parameters - these are set up by R in the function 'writeReport' when copying the 
//   template for this script from arrayQualityMetrics/inst/scripts into the report.

var highlightInitial = ["false","false","true","false","false","false","true","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","true","false","false","false","false","false","true","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false","false"];
var arrayMetadata    = [{"array":"1","sampleNames":"GSM1256653_123.CEL.gz","index":"1"},{"array":"2","sampleNames":"GSM1256654_118.CEL.gz","index":"2"},{"array":"3","sampleNames":"GSM1256655_103.CEL.gz","index":"3"},{"array":"4","sampleNames":"GSM1256656_132.CEL.gz","index":"4"},{"array":"5","sampleNames":"GSM1256657_119.CEL.gz","index":"5"},{"array":"6","sampleNames":"GSM1256658_133.CEL.gz","index":"6"},{"array":"7","sampleNames":"GSM1256659_134.CEL.gz","index":"7"},{"array":"8","sampleNames":"GSM1256660_104.CEL.gz","index":"8"},{"array":"9","sampleNames":"GSM1256661_121.CEL.gz","index":"9"},{"array":"10","sampleNames":"GSM1256662_135.CEL.gz","index":"10"},{"array":"11","sampleNames":"GSM1256663_108.CEL.gz","index":"11"},{"array":"12","sampleNames":"GSM1256664_109.CEL.gz","index":"12"},{"array":"13","sampleNames":"GSM1256665_136.CEL.gz","index":"13"},{"array":"14","sampleNames":"GSM1256666_102.CEL.gz","index":"14"},{"array":"15","sampleNames":"GSM1256667_122.CEL.gz","index":"15"},{"array":"16","sampleNames":"GSM1256668_137.CEL.gz","index":"16"},{"array":"17","sampleNames":"GSM1256669_138.CEL.gz","index":"17"},{"array":"18","sampleNames":"GSM1256670_139.CEL.gz","index":"18"},{"array":"19","sampleNames":"GSM1256671_140.CEL.gz","index":"19"},{"array":"20","sampleNames":"GSM1256672_94.CEL.gz","index":"20"},{"array":"21","sampleNames":"GSM1256673_110.CEL.gz","index":"21"},{"array":"22","sampleNames":"GSM1256674_95.CEL.gz","index":"22"},{"array":"23","sampleNames":"GSM1256675_90.CEL.gz","index":"23"},{"array":"24","sampleNames":"GSM1256676_141.CEL.gz","index":"24"},{"array":"25","sampleNames":"GSM1256677_89.CEL.gz","index":"25"},{"array":"26","sampleNames":"GSM1256678_105.CEL.gz","index":"26"},{"array":"27","sampleNames":"GSM1256679_106.CEL.gz","index":"27"},{"array":"28","sampleNames":"GSM1256680_107.CEL.gz","index":"28"},{"array":"29","sampleNames":"GSM1256681_120.CEL.gz","index":"29"},{"array":"30","sampleNames":"GSM1256682_78.CEL.gz","index":"30"},{"array":"31","sampleNames":"GSM1256683_100.CEL.gz","index":"31"},{"array":"32","sampleNames":"GSM1256684_79.CEL.gz","index":"32"},{"array":"33","sampleNames":"GSM1256685_72.CEL.gz","index":"33"},{"array":"34","sampleNames":"GSM1256686_82.CEL.gz","index":"34"},{"array":"35","sampleNames":"GSM1256687_124.CEL.gz","index":"35"},{"array":"36","sampleNames":"GSM1256688_91.CEL.gz","index":"36"},{"array":"37","sampleNames":"GSM1256689_80.CEL.gz","index":"37"},{"array":"38","sampleNames":"GSM1256690_96.CEL.gz","index":"38"},{"array":"39","sampleNames":"GSM1256691_73.CEL.gz","index":"39"},{"array":"40","sampleNames":"GSM1256692_92.CEL.gz","index":"40"},{"array":"41","sampleNames":"GSM1256693_74.CEL.gz","index":"41"},{"array":"42","sampleNames":"GSM1256694_97.CEL.gz","index":"42"},{"array":"43","sampleNames":"GSM1256695_125.CEL.gz","index":"43"},{"array":"44","sampleNames":"GSM1256696_130.CEL.gz","index":"44"},{"array":"45","sampleNames":"GSM1256697_142.CEL.gz","index":"45"},{"array":"46","sampleNames":"GSM1256698_81.CEL.gz","index":"46"},{"array":"47","sampleNames":"GSM1256699_111.CEL.gz","index":"47"},{"array":"48","sampleNames":"GSM1256700_98.CEL.gz","index":"48"},{"array":"49","sampleNames":"GSM1256701_126.CEL.gz","index":"49"},{"array":"50","sampleNames":"GSM1256702_75.CEL.gz","index":"50"},{"array":"51","sampleNames":"GSM1256703_113.CEL.gz","index":"51"},{"array":"52","sampleNames":"GSM1256704_143.CEL.gz","index":"52"},{"array":"53","sampleNames":"GSM1256705_85.CEL.gz","index":"53"},{"array":"54","sampleNames":"GSM1256706_76.CEL.gz","index":"54"},{"array":"55","sampleNames":"GSM1256707_114.CEL.gz","index":"55"},{"array":"56","sampleNames":"GSM1256708_88.CEL.gz","index":"56"},{"array":"57","sampleNames":"GSM1256709_112.CEL.gz","index":"57"},{"array":"58","sampleNames":"GSM1256710_131.CEL.gz","index":"58"},{"array":"59","sampleNames":"GSM1256711_144.CEL.gz","index":"59"},{"array":"60","sampleNames":"GSM1256712_127.CEL.gz","index":"60"},{"array":"61","sampleNames":"GSM1256713_83.CEL.gz","index":"61"},{"array":"62","sampleNames":"GSM1256714_145.CEL.gz","index":"62"},{"array":"63","sampleNames":"GSM1256715_93.CEL.gz","index":"63"},{"array":"64","sampleNames":"GSM1256716_128.CEL.gz","index":"64"},{"array":"65","sampleNames":"GSM1256717_129.CEL.gz","index":"65"},{"array":"66","sampleNames":"GSM1256718_77.CEL.gz","index":"66"},{"array":"67","sampleNames":"GSM1256719_115.CEL.gz","index":"67"},{"array":"68","sampleNames":"GSM1256720_1.CEL.gz","index":"68"},{"array":"69","sampleNames":"GSM1256721_21.CEL.gz","index":"69"},{"array":"70","sampleNames":"GSM1256722_2.CEL.gz","index":"70"},{"array":"71","sampleNames":"GSM1256723_3.CEL.gz","index":"71"},{"array":"72","sampleNames":"GSM1256724_22.CEL.gz","index":"72"},{"array":"73","sampleNames":"GSM1256725_4.CEL.gz","index":"73"},{"array":"74","sampleNames":"GSM1256726_23.CEL.gz","index":"74"},{"array":"75","sampleNames":"GSM1256727_24.CEL.gz","index":"75"},{"array":"76","sampleNames":"GSM1256728_27.CEL.gz","index":"76"},{"array":"77","sampleNames":"GSM1256729_5.CEL.gz","index":"77"},{"array":"78","sampleNames":"GSM1256730_25.CEL.gz","index":"78"},{"array":"79","sampleNames":"GSM1256731_6.CEL.gz","index":"79"},{"array":"80","sampleNames":"GSM1256732_28.CEL.gz","index":"80"},{"array":"81","sampleNames":"GSM1256733_7.CEL.gz","index":"81"},{"array":"82","sampleNames":"GSM1256734_8.CEL.gz","index":"82"},{"array":"83","sampleNames":"GSM1256735_56.CEL.gz","index":"83"},{"array":"84","sampleNames":"GSM1256736_57.CEL.gz","index":"84"},{"array":"85","sampleNames":"GSM1256737_58.CEL.gz","index":"85"},{"array":"86","sampleNames":"GSM1256738_59.CEL.gz","index":"86"},{"array":"87","sampleNames":"GSM1256739_35.CEL.gz","index":"87"},{"array":"88","sampleNames":"GSM1256740_60.CEL.gz","index":"88"},{"array":"89","sampleNames":"GSM1256741_36.CEL.gz","index":"89"},{"array":"90","sampleNames":"GSM1256742_37.CEL.gz","index":"90"},{"array":"91","sampleNames":"GSM1256743_61.CEL.gz","index":"91"},{"array":"92","sampleNames":"GSM1256744_38.CEL.gz","index":"92"},{"array":"93","sampleNames":"GSM1256745_39.CEL.gz","index":"93"},{"array":"94","sampleNames":"GSM1256746_62.CEL.gz","index":"94"},{"array":"95","sampleNames":"GSM1256747_63.CEL.gz","index":"95"},{"array":"96","sampleNames":"GSM1256748_40.CEL.gz","index":"96"},{"array":"97","sampleNames":"GSM1256749_64.CEL.gz","index":"97"},{"array":"98","sampleNames":"GSM1256750_65.CEL.gz","index":"98"},{"array":"99","sampleNames":"GSM1256751_50.CEL.gz","index":"99"},{"array":"100","sampleNames":"GSM1256752_51.CEL.gz","index":"100"},{"array":"101","sampleNames":"GSM1256753_41.CEL.gz","index":"101"},{"array":"102","sampleNames":"GSM1256754_42.CEL.gz","index":"102"},{"array":"103","sampleNames":"GSM1256755_66.CEL.gz","index":"103"},{"array":"104","sampleNames":"GSM1256756_43.CEL.gz","index":"104"},{"array":"105","sampleNames":"GSM1256757_67.CEL.gz","index":"105"},{"array":"106","sampleNames":"GSM1256758_68.CEL.gz","index":"106"},{"array":"107","sampleNames":"GSM1256759_44.CEL.gz","index":"107"},{"array":"108","sampleNames":"GSM1256760_45.CEL.gz","index":"108"},{"array":"109","sampleNames":"GSM1256761_46.CEL.gz","index":"109"},{"array":"110","sampleNames":"GSM1256762_47.CEL.gz","index":"110"},{"array":"111","sampleNames":"GSM1256763_48.CEL.gz","index":"111"},{"array":"112","sampleNames":"GSM1256764_52.CEL.gz","index":"112"},{"array":"113","sampleNames":"GSM1256765_53.CEL.gz","index":"113"},{"array":"114","sampleNames":"GSM1256766_29.CEL.gz","index":"114"},{"array":"115","sampleNames":"GSM1256767_30.CEL.gz","index":"115"},{"array":"116","sampleNames":"GSM1256768_31.CEL.gz","index":"116"},{"array":"117","sampleNames":"GSM1256769_32.CEL.gz","index":"117"},{"array":"118","sampleNames":"GSM1256770_9.CEL.gz","index":"118"},{"array":"119","sampleNames":"GSM1256771_10.CEL.gz","index":"119"},{"array":"120","sampleNames":"GSM1256772_11.CEL.gz","index":"120"},{"array":"121","sampleNames":"GSM1256773_148.CEL.gz","index":"121"},{"array":"122","sampleNames":"GSM1256774_86.CEL.gz","index":"122"},{"array":"123","sampleNames":"GSM1256775_99.CEL.gz","index":"123"},{"array":"124","sampleNames":"GSM1256776_87.CEL.gz","index":"124"},{"array":"125","sampleNames":"GSM1256777_116.CEL.gz","index":"125"},{"array":"126","sampleNames":"GSM1256778_84.CEL.gz","index":"126"},{"array":"127","sampleNames":"GSM1256779_101.CEL.gz","index":"127"},{"array":"128","sampleNames":"GSM1256780_146.CEL.gz","index":"128"},{"array":"129","sampleNames":"GSM1256781_147.CEL.gz","index":"129"},{"array":"130","sampleNames":"GSM1256782_117.CEL.gz","index":"130"},{"array":"131","sampleNames":"GSM1256783_70.CEL.gz","index":"131"},{"array":"132","sampleNames":"GSM1256784_12.CEL.gz","index":"132"},{"array":"133","sampleNames":"GSM1256785_13.CEL.gz","index":"133"},{"array":"134","sampleNames":"GSM1256786_14.CEL.gz","index":"134"},{"array":"135","sampleNames":"GSM1256787_26.CEL.gz","index":"135"},{"array":"136","sampleNames":"GSM1256788_15.CEL.gz","index":"136"},{"array":"137","sampleNames":"GSM1256789_69.CEL.gz","index":"137"},{"array":"138","sampleNames":"GSM1256790_54.CEL.gz","index":"138"},{"array":"139","sampleNames":"GSM1256791_33.CEL.gz","index":"139"},{"array":"140","sampleNames":"GSM1256792_55.CEL.gz","index":"140"},{"array":"141","sampleNames":"GSM1256793_16.CEL.gz","index":"141"},{"array":"142","sampleNames":"GSM1256794_49.CEL.gz","index":"142"},{"array":"143","sampleNames":"GSM1256795_71.CEL.gz","index":"143"},{"array":"144","sampleNames":"GSM1256796_17.CEL.gz","index":"144"},{"array":"145","sampleNames":"GSM1256797_34.CEL.gz","index":"145"},{"array":"146","sampleNames":"GSM1256798_18.CEL.gz","index":"146"},{"array":"147","sampleNames":"GSM1256799_19.CEL.gz","index":"147"},{"array":"148","sampleNames":"GSM1256800_20.CEL.gz","index":"148"}];
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
