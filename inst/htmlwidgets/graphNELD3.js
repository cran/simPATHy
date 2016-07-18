HTMLWidgets.widget({

  name: 'graphNELD3',

  type: 'output',

  initialize: function(el, width, height) {

    d3.select(el).append("svg")
        .attr("width", width)
        .attr("height", height);

    return d3.layout.force();
  },
  resize: function(el, width, height, force) {


    d3.select(el).select("svg")
        .attr("width", width)
        .attr("height", height);

    force.size([width, height]).resume();
  },

  renderValue: function(el, x, force) {

    var drag = force.drag()
        .on("dragstart", dragstart);

    function dragstart(d) {
        d3.event.sourceEvent.preventDefault();
        d3.event.sourceEvent.stopPropagation();
    }

    var svg = d3.select(el).select("svg");
              svg.selectAll("*").remove();

    svg = d3.select(el).select("svg")
              .append("g").attr("class","zoom-layer");

    var zoom = d3.behavior.zoom().on("zoom",redraw);
    function redraw() {
        d3.select(el).select(".zoom-layer").attr("transform",
              "translate(" + d3.event.translate + ")"+
              " scale(" + d3.event.scale + ")");
    }

    d3.select(el).select("svg")
      .attr("pointer-events", "all")
      .call(zoom);

		var dataNode=HTMLWidgets.dataframeToD3(x.nodes);
		var dataEdge=HTMLWidgets.dataframeToD3(x.edges);

    // opzioni in ingresso
    var type = x.options.type;
    var colNode = x.options.colNode;
    var legend = x.options.legend;

    var width = el.offsetWidth;
    var height = el.offsetHeight;

    var colors = d3.scale.category10();

    // create d3 force layout
    force
      .nodes(dataNode)
      .links(dataEdge)
      .size([width, height])
      .linkDistance(50)
      .charge(-300)
      .start();

      var border = svg.selectAll("border")
            .data(dataEdge)
            .attr("id",function(d,i){ return i+100 })
				    .enter()
			    	.append("line")
				    .style("stroke","#ccc")
				    .style("stroke-width", 2)
				    .style("opacity",function(d){
				        if(d.colour=="#c0c0c0") return 0.2;
				        else return 1;
				    });

			var edges = svg.selectAll("link")
				.data(dataEdge)
				.attr("id",function(d,i){ return i })
				.enter()
				.append("line")
				.style("stroke",function(d){ return d.colour;})
				.style("stroke-width", 1.5)
				.style("stroke-dasharray",function(d){
				  if(d.colour=="#c0c0c0") return "1,2";
				  else return "null";
				});

      if(x.options.edgeClick){
        edges.on("click",edgeClick);
      }

			if( type == "dag" ){
			  edges.style("marker-end",function(d,i){ return "url(#"+d.source.name+"_"+d.target.name+")"  });
			}

      var nodes = svg.selectAll(".node")
            .data(dataNode)
            .enter().append("g")
            .attr("class","node")
            .style("opacity",1)
            .call(force.drag);

      nodes.append("circle")
        .attr("r", 10)
        .style("fill",colNode)
        .style("stroke", "white");

      nodes.append("svg:text")
        .attr("class", "nodetext")
        .attr("dx", 13)
        .attr("dy", ".35em")
        .text(function(d) { return d.name })
        .style("font", "5px")
        .style("fill","#6e6e6e");

			var defs = svg.append("defs");


			var M = defs.selectAll("marker");
      M.data(dataEdge)
			 .enter().append("marker")
			 .attr("id",function(d,i){ return d.source.name+"_"+d.target.name })
			 .attr("viewBox","0 -5 10 10")
			 .attr("refX",20.5)
			 .attr("refY",0)
			 .attr("markerWidth",6)
			 .attr("markerHeight",6)
			 .attr("orient","auto")
			 .append("path")
			 .attr("d", "M0,-5L10,0L0,5")
       .style("fill",function(d,i){ return d.colour })
       .style("stroke","#ccc")
       .style("stroke-width",0.2);


			M.data(["clicked"])
			 .enter().append("marker")
			 .attr("id",function(d){ return d })
			 .attr("viewBox","0 -5 10 10")
			 .attr("refX",20)
			 .attr("refY",0)
			 .attr("markerWidth",4)
			 .attr("markerHeight",4)
			 .attr("orient","auto")
			 .append("path")
			 .attr("d", "M0,-5L10,0L0,5")
       .style("fill","black")
       .style("stroke", "black");

			var SS=[];
			var TT=[];

      function edgeClick(d){

        var s=d.source.index+1;
        var t=d.target.index+1;
        var conta=0;
        var index=[];

        if(SS.length > 0){
          for(var i=0; i< SS.length; i++){
            if( SS[i]==s & TT[i]==t ) index.push(i);
          }
        }

        if( index.length === 0 ){
          // lo evidenzio e lo aggiungo a SS e TT
          d3.select(this)
            .style("stroke", "black")
            .style("stroke-width", 2.5)
            .style("marker-end",function(d) {
              if(type=="dag")
                return "url(#clicked)";
            });


          SS.push(s);
          TT.push(t);

        } else {

          d3.select(this)
            .style("stroke", function(d){ return d.colour; })
			    	.style("stroke-width", 1.5)
			    	.style("marker-end",function(d) {
              if(type=="dag")
                return  "url(#"+d.source.name+"_"+d.target.name+")" ;//"url(#end)";
            });

          SS.splice(index,1);
          TT.splice(index,1);
        }

        nodes.selectAll("circle").style("fill",function(n){
            var presS=SS.indexOf(n.index+1);
            var presT=TT.indexOf(n.index+1);

            if( presS=== -1 & presT=== -1 ) return colNode;
            else return "#6495ed";

        });

        var EE=[];
        EE.push(SS);
        EE.push(TT);

        Shiny.onInputChange("mydata",EE);

      }


			force.on("tick", function() {

				edges.attr("x1", function(d) { return d.source.x; })
					 .attr("y1", function(d) { return d.source.y; })
					 .attr("x2", function(d) { return d.target.x; })
					 .attr("y2", function(d) { return d.target.y; });

        border.attr("x1", function(d) { return d.source.x; })
					 .attr("y1", function(d) { return d.source.y; })
					 .attr("x2", function(d) { return d.target.x; })
					 .attr("y2", function(d) { return d.target.y; });

			  nodes.attr("transform",function(d){
			      return "translate(" + d.x + "," + d.y + ")";
			  });


			});


      if(legend){

          var linearGradient = defs.append("linearGradient")
                 .attr("id","linear-gradient");

          linearGradient.append("stop")
               .attr("offset","0%")
               .attr("stop-color","#8B0000");
          linearGradient.append("stop")
               .attr("offset","50%")
               .attr("stop-color","white");
          linearGradient.append("stop")
               .attr("offset","100%")
               .attr("stop-color","#000080");

          var recLegend = d3.select(el).select("svg")
                            .append("g").attr("class","legend");

          //var wleg=300;
          var wleg= 300;
          var start=(width-wleg)/2;
          var step = wleg/(x.dataLegend.value.length-1);

          // add rect
          recLegend.append("rect")
             .attr("width",wleg)
             .attr("height",10)
             .attr("x",start )
             .attr("y", height-20 )
             .style("fill","url(#linear-gradient)");

          recLegend.selectAll("text-legend")
              .data(x.dataLegend.value).enter()
              .append("text")
              .attr("x", function(d,i){ return start+step*i } )
              .attr("y", height -25)
              .text(function(d,i){ return d })
              .style("text-anchor","middle")
              .style("font-size", "10px")
              .style("fill","#6e6e6e");


      }



  }


});
