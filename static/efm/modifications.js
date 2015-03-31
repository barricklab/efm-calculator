// add parser through the tablesorter addParser method
// Allows us to sort rates even if they are in scientific notation
$.tablesorter.addParser({
    // set a unique id
    id: 'scinot',
    is: function(s) {
        // Matches scientific notation including leading "~", "<", or ">" chars
        return /[+\-\>\<\~]?[\s]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)?/.test(s);
    },
    format: function(s) {
        // Strip out leading characters and convert dashes to zero
        s = s.replace("<","").replace("~","").replace(/^[\u2013|\u2014|-]$/, "0");
        return $.tablesorter.formatFloat(s);
    },
    type: 'numeric'
});

$(function() {
    $("#mytable").tablesorter({
        headers : {
            5 : { sorter: 'scinot' }},
        sortList: [[5,1]]
    });
});

function getRandColor() {
    var r = function () { return Math.floor(Math.random()*256) };
    return "rgb(" + r() + "," + r() + "," + r() + ")";
}

function getRateColor(rate) {
    var pct = (-Math.log(rate)/Math.LN10)/6;
    var percentColors = [
        { pct: 0.6, color: { r: 220, g: 20, b: 60 } },
        { pct: 1.0, color: { r: 255, g: 255, b: 0 } } ];

    for (var i = 1; i < percentColors.length - 1; i++) {
        if (pct < percentColors[i].pct) {
            break;
        }
    }
    var lower = percentColors[i - 1];
    var upper = percentColors[i];
    var range = upper.pct - lower.pct;
    var rangePct = (pct - lower.pct) / range;
    var pctLower = 1 - rangePct;
    var pctUpper = rangePct;
    var color = {
        r: Math.floor(lower.color.r * pctLower + upper.color.r * pctUpper),
        g: Math.floor(lower.color.g * pctLower + upper.color.g * pctUpper),
        b: Math.floor(lower.color.b * pctLower + upper.color.b * pctUpper)
    };
    return 'rgb(' + [color.r, color.g, color.b].join(',') + ')';
}

// This function sets a minimum glyph size so that nothing disappears when you zoom out

function setMinGlyphSize(feature, length) {
    if((feature.lane.chart.pixelsToNts() * feature.length) < length) {
          // create rectangle with the same characteristics and attach it to the same chart
          var r = new Rect(feature.type, feature.position, feature.lane.chart.ntsToPixels() * length, feature.opts);
          r.lane = feature.lane;
          r.ctx = feature.ctx;
          r._draw();
          // return true to stop normal drawing of glyph
          return true;
    }
    // return false to allow normal drawing
    return false;
}


var custom_Scribl = Scribl.extend({initScrollable: function() {
      var scrollStartMin;

      if (!this.scrolled){
         // create divs
         var parentDiv = document.createElement('div');
         var canvasContainer = document.createElement('div');
         var sliderDiv = document.createElement('div');
         sliderDiv.id = 'scribl-zoom-slider';
         sliderDiv.className = 'slider';
         sliderDiv.style.cssFloat = 'left';
         sliderDiv.style.height = (new String(this.canvas.height * .75)) + 'px';
         sliderDiv.style.margin = '30px auto auto -20px'

         // grab css styling from canavs
         parentDiv.style.cssText = this.canvas.style.cssText;
         this.canvas.style.cssText = '';
         parentWidth = parseInt(this.canvas.width) + 25;
         parentDiv.style.width = parentWidth + 'px';
         canvasContainer.style.width = this.canvas.width + 'px';
         canvasContainer.style.overflow = 'auto';
         canvasContainer.id = 'scroll-wrapper';



         this.canvas.parentNode.replaceChild(parentDiv, this.canvas);
         parentDiv.appendChild(sliderDiv);
         canvasContainer.appendChild(this.canvas);
         parentDiv.appendChild(canvasContainer);

         jQuery(canvasContainer).dragscrollable({dragSelector: 'canvas:first', acceptPropagatedEvent: false});
      }

      var totalNts =  this.scale.max - this.scale.min;
      var scrollStartMax = this.scrollValues[1] || this.scale.max;
      if( this.scrollValues[0] != undefined)
          scrollStartMin = this.scrollValues[0];
      else
          scrollStartMin = this.scale.min;

      var viewNts = scrollStartMax - scrollStartMin;
      var viewNtsPerPixel = viewNts / document.getElementById('scroll-wrapper').style.width.split('px')[0];

      var canvasWidth = (totalNts / viewNtsPerPixel) || 100;
      this.canvas.width = canvasWidth;
      this.width = canvasWidth - 30;
      schart = this;
      var zoomValue = (scrollStartMax - scrollStartMin) / (this.scale.max - this.scale.min) * 100 || 1;
      var minValue = (document.getElementById('scroll-wrapper').style.width.split('px')[0]/totalNts)*6;

      jQuery(sliderDiv).slider({
         orientation: 'vertical',
         range: 'min',
         min: minValue,
         max: 100,
         value: zoomValue,
         slide: function( event, ui ) {
            var totalNts = schart.scale.max - schart.scale.min;
            var width = ui.value / 100 * totalNts;
            var widthPixels = ui.value / 100 * schart.canvas.width;
            var canvasContainer = document.getElementById('scroll-wrapper');
            var center = canvasContainer.scrollLeft + parseInt(canvasContainer.style.width.split('px')[0]) / 2;

            // get min max pixels
            var minPixel = center - widthPixels/2;
            var maxPixel = center + widthPixels/2;

            // convert to nt
            var min = schart.scale.min + (minPixel / schart.canvas.width) * totalNts;
            var max = schart.scale.min + (maxPixel / schart.canvas.width) * totalNts;

            schart.scrollValues = [min, max];
            schart.ctx.clearRect(0, 0, schart.canvas.width, schart.canvas.height);
            schart.draw();
         }
      });


      var startingPixel = (scrollStartMin - this.scale.min) / totalNts * this.canvas.width;
      document.getElementById('scroll-wrapper').scrollLeft = startingPixel;
      this.scrolled = true;
    }
});