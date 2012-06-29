jQuery(document).ready(function(){
 $("#hidePanel").click(function(){
 $("#panel").animate({marginLeft:"-140px"}, 500 );
 $("#colleft").animate({width:"0px", opacity:0}, 400 );
 $("#showPanel").show("normal").animate({width:"28px", opacity:1}, 300);
 //$("#colright").animate({marginLeft:"50px"}, 500);
 });
 $("#showPanel").click(function(){
 //$("#colright").animate({marginLeft:"200px"}, 200);
 $("#panel").animate({marginLeft:"0px"}, 400 );
 $("#colleft").animate({width:"140px", opacity:1}, 400 );
 $("#showPanel").animate({width:"0px", opacity:0}, 600).hide("slow");
 });
});
