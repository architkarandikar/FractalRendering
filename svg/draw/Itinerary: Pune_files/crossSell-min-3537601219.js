var CROSSSELL=(function(h){var v=h("#overlay_rcc");var b=h("#pop_up_modal_rcc");var e=h(".timer_rcc").length!==0;var x="Car.XSell.V3";var A=null;var t=null;var j=false;var m=false;var c=0;function f(F,E,C,D){h.ajax({type:"POST",url:"https://"+window.location.host+"/Confirmation-Flight/carReservationCrossSellReserve?tripid="+E,data:F,mimeType:"application/json",success:function(G){if(G){for(var H=0;H<=2;H++){if(D===C[H]){s_exp.pageName="page.Flight.Interstitial.Car.XSell.Confirmation";s_exp.prop16=x;s_exp.eVar28=x;s_exp.t();h("#loading_rcc").remove();h("#one_click_confirmation_container_rcc").show();h("#one_click_confirmation_rcc_"+(H+1)).show();}}}else{h("#loading_rcc").remove();h("#one_click_error_rcc").show();s_exp.pageName="page.Flight.Interstitial.Car.Error";s_exp.prop16=x;s_exp.eVar28=x;s_exp.t();}},error:function(G,I,H){h("#loading_rcc").remove();h("#one_click_error_rcc").show();s_exp.pageName="page.Flight.Interstitial.Car.Error";s_exp.prop16=x;s_exp.eVar28=x;s_exp.t();}});}function i(){var C=h("html").hasClass("ie8");if(A==null){A=function(G){var I=h("#in_page_checkout_div");if(b.hasClass("inpagecko")){if(G.data){try{var H=JSON.parse(G.data);if(H["action"]==="changeheight"){if(H["height"]){var D=H["height"];var J=(h(window).width()-900)/2;if(J<0){J=20;}var F="75em";if(C){F="80em";}b.css({left:J,top:"0px",padding:"0px","height":D,width:F});I.css({visibility:"visible"});}}else{if(H["action"]==="changeboth"){if(H["height"]&&H["width"]){var F=H["width"];var D=H["height"];var J=(h(window).width()-F)/2;if(J<0){J=20;}b.css({left:J,top:"0px",padding:"0px","height":D,width:F});I.css({visibility:"visible"});}}else{if(H["action"]==="cancel"){p();}else{if(H["action"]==="hidecontent"){I.css({visibility:"hidden"});}else{if(H["action"]==="showcontent"){I.css({visibility:"visible"});}else{if(H["action"]==="blockhide"){j=true;}else{if(H["action"]==="unblockhide"){j=false;}else{if(H["action"]==="complete"){m=true;}}}}}}}}}catch(E){}}}};if(window.addEventListener){window.addEventListener("message",A);}else{window.attachEvent("onmessage",A);}h("#in_page_checkout").load(function(){c++;var D=c;window.setTimeout(function(){if(D!=c){return;}var I=h("#in_page_checkout_div");if(b.hasClass("inpagecko")){if(I.is(":visible")===false){var E=400;var G=900;try{var H=h("#in_page_checkout").contents();E=H.find("body").height();G=H.find("div.site-header").width();}catch(F){}window.postMessage(JSON.stringify({action:"changeboth",width:G,height:E}),"*");}if(j){window.postMessage(JSON.stringify({action:"unblockhide"}),"*");}}},2000);});}}function p(){var D,C=h("#in_page_checkout_div");if(j===false){C.css({visibility:"hidden"});b.removeClass("inpagecko");D=(h(window).width()-t.width)/2;b.css({left:D,top:t.top,padding:t.padding,width:t.width,height:t.height});h("#one_click_offer_rcc").fadeIn();}}function g(G){var F,D=h("#in_page_checkout_div"),E=h("#in_page_checkout"),C="https://"+window.location.host+G.data("path")+"car.id="+G.data("piid")+"&previousTripId="+G.data("tripid")+"&inPageCko=true"+"&totalPriceShown="+G.data("totalpriceshown");if(t===null){t=h.extend(t,b.position(),{width:b.width()},{height:b.height()},{padding:b.css("padding")});}D.css({visibility:"hidden"});E.attr("src","");h("#one_click_offer_rcc").fadeOut({complete:function(){b.addClass("inpagecko");b.css({height:t.height+36,width:t.width+36,padding:"0px"});i();E.attr("src",C);}});}function n(){var E=h("#oneClickReservation").html();var H=h("#seamlessCheckoutApi").html();var I=h("#seamlessCKOInPageCKO").html();b.css("left",(h(window).width()-b.width())/2+"px");h("#reservation_cross_sell_more_car_rcc").on("click",function(){s_exp.clearOmnitureObject();s_exp.prop16="CKO.FLT.XSell.SeeMoreCars";s_exp.eVar28="CKO.FLT.XSell.SeeMoreCars";s_exp.t();D();});h(".skip_offer_close_x_rcc, .close_offer_rcc").on("click",function(){s_exp.clearOmnitureObject();s_exp.prop16="Car.XSell.Close";s_exp.eVar28="Car.XSell.Close";s_exp.t();D();});v.bind("click",function(){D();});var D=function(){var L=h("#loading_rcc");if(m===true){v.remove();b.remove();}else{if(b.hasClass("inpagecko")){p();}else{if(L.length===0||!L.is(":visible")){v.remove();b.remove();}}}};var K=function(){s_exp.clearOmnitureObject();s_exp.prop16="CKO.FLT.XSell.CarSelected";s_exp.eVar28="CKO.FLT.XSell.CarSelected";s_exp.t();};h(".button_rcc").on("click",function(){K();if(E=="true"){var L=this,Q=h(L);if(H=="true"&&Q.data("merchant")){if(I==="1"){d("7280","1");g(Q);}else{d("7280","0");u(Q);}}else{if(e){clearInterval(F);}var P=null;if(Q.data("merchant")){P=window.open();}h("#one_click_offer_rcc").hide();Crosssell.completeBooking(Q.data("tripid"),Q.data("piid"),Q.data("path"),"car",Q.data("merchant"),h("#loading_rcc"),function(S){if(S.hasOwnProperty("redirectUrl")){q(P,S);}else{o(L,S);}},function(T){if(P!=null){P.close();}h("#one_click_error_rcc").show();s_exp.clearOmnitureObject();if(T.omnitureData!==undefined){var S=h.parseJSON(T.omnitureData);h.each(S,function(U,V){if(h.isPlainObject(V)){h.each(V,function(W,X){s_exp[W]=X;});}else{s_exp[U]=V;}});}s_exp.pageName="page.Flight.Interstitial.Car.Error";s_exp.prop16=x;s_exp.eVar28=x;window.s_exp_ado=window.s_exp_ado||[];window.s_exp_ado=[["TRL",s_exp.prop71],["OrderNumber",s_exp.prop72]];s_exp.t();});}}else{if(e){clearInterval(F);}var O=h(this).parents(".car_info_box_rcc");var N=O.find(".trip_id_rcc").text();var M=h(".button_rcc");var R={carClass:O.find(".car_category_rcc").text(),quotedPrice:O.find(".actual_price_rcc").text().replace(/[$,]/g,""),airportCode:O.find(".airport_code_rcc").text(),carModel:O.find(".model_rcc").text(),vendorName:O.find(".vendor_name_rcc").text(),pickUp:O.find(".pick_up_rcc").text().split(" ").join(""),dropOff:O.find(".drop_off_rcc").text().split(" ").join("")};h("#one_click_offer_rcc").hide();h("#loading_rcc").show();f(R,N,M,this);}});v.show();h("#pop_up_modal_rcc").show();s_exp.clearOmnitureObject();s_exp.pageName="page.Flight.Interstitial.Car.XSell";s_exp.prop16=x;s_exp.eVar28=x;s_exp.t();if(e){var F;var J;function C(){J=h(".time_minutes_rcc").html()*60*1000;F=self.setInterval(function(){G();},1000);}function G(){J=J-1000;var L=Math.floor(J/60000);var M=Math.floor(J%60000/1000);if(M<10){M="0"+M;}if(J===0){clearInterval(F);h("#pop_up_modal_rcc").remove();h("#overlay_rcc").remove();}else{h(".time_minutes_rcc").html(L);h(".time_seconds_rcc").html(M);}}C();}}function y(E,C){var D=E+"."+C;if(typeof s_exp.prop34==="undefined"||s_exp.prop34===""){s_exp.prop34=D;}else{if(s_exp.prop34.indexOf(D)==(-1)){s_exp.prop34=s_exp.prop34+"|"+D;}}s_exp.tl(this,"o","ape:CKO:FlightConfirmation.InPageCKO");}function d(F,E){y(F,E);try{var D='{"guid":"'+this.guid+'", "tpid":'+tpid+', "eapid":'+eapid+' , "evaluatedExperiments" : [{ "id": '+F+', "value": '+E+'}], "pageName":"page.Flight.Checkout.Confirmation"}';h.ajax({url:"/api/bucketing/v1/logExperiments",dataType:"json",data:D,type:"POST",contentType:"application/json",cache:false,success:function(G){if(G){if(G.successes!==null&&G.successes.length>0){h("#in_page_checkout").attr("data-logresult",G.successes[0].id);}}},error:function(){setTimeout(function(){throw"LogExperiment for 7280 failed";},100);}});}catch(C){setTimeout(function(){throw"LogExperiment for 7280 failed";},100);}}function u(D){var C=D.data("path");if(C.indexOf("c=")>0){var E=C.substring(C.indexOf("=")+1,C.indexOf("&"));h("#c").val(E);}h("#carId").val(D.data("piid"));h("#previousTripId").val(D.data("tripid"));h("#fallbackUrl").val(encodeURIComponent(D.data("fallbackurl")));h("#totalPriceShown").val(D.data("totalpriceshown"));h("#seamless_api_form").attr("action",C.substring(0,C.indexOf("?")));h("#seamless_api_form").submit();}var o=function(E,G){var C=h(".button_rcc");for(var D=0;D<=2;D++){if(E===C[D]){h("#one_click_confirmation_container_rcc").show();h("#one_click_confirmation_rcc_"+(D+1)).show();h("#itinerary_link_"+(D+1)).attr("href",G.itineraryUrl);h("#itinerary_link_"+(D+1)).text(h("#itinerary_link_"+(D+1)).text()+G.itineraryNumber);}}if(G.omnitureData!==undefined){var F=h.parseJSON(G.omnitureData);s_exp.clearOmnitureObject();
h.each(F,function(H,I){if(h.isPlainObject(I)){h.each(I,function(J,K){s_exp[J]=K;});}else{s_exp[H]=I;}});}s_exp.pageName="page.Flight.Interstitial.Car.XSell.Confirmation";s_exp.prop16=x;s_exp.eVar28=x;window.s_exp_ado=window.s_exp_ado||[];window.s_exp_ado=[["TRL",s_exp.prop71],["OrderNumber",s_exp.prop72]];s_exp.t();};var q=function(D,C){h("#one_click_offer_rcc").show();D.location.replace(C.redirectUrl);};if(h("#pop_up_modal_rcc").length!==0){var a=0;if(typeof flightConfirmationCrossSellDelay!=="undefined"){a=flightConfirmationCrossSellDelay;}setTimeout(function(){n();},a);}var B="CROSSSELL.getRecommendations()::";var w=function(C){return decodeURI((RegExp(C+"="+"(.+?)(&|$)").exec(location.search)||[,null])[1]);};var s=function(C){if(C){C=C.replace(/http:\/\//,"https://").replace(/l.jpg/i,"t.jpg");}else{C="https://"+window.location.host+"/static/default/default/images/hotResult/noPhotosAvailSmall.gif";}return C;};var z=function(E){var C=Math.round(E*10)/10;var D=String(C);if(D.length>1){D=D.replace(/\./,"-");}else{D+="-0";}return"value-"+D;};var r=function(C){C=C.replace("https","http");return C;};var l=function(D,M,J){if(D&&D.length){var N,C,H,F,I,K,E,L=(D.length>=3)?3:D.length;for(var G=0;G<L;G++){N=h("#"+M).clone();N.removeClass("hidden");H=D[G];F=N.find("span.stars0-0");F.removeClass("stars0-0");F.addClass("stars"+H.starRating.replace(/\./,"-"));F.html(H.starRating);K=(H.reviewTotal>0)?true:false;if(K){N.find("#hotel_cross_sell_box_review").removeClass("hidden");I=N.find("span.value-0-0");I.removeClass("value-0-0");I.addClass(z(H.reviewRating));E=N.find(".handleSingularAndPlural").html().split("|");if(H.reviewTotal==1){N.find(".handleSingularAndPlural").html(E[0].replace(/^\s\s*/,"").replace(/\s\s*$/,""));}else{N.find(".handleSingularAndPlural").html(E[1].replace(/^\s\s*/,"").replace(/\s\s*$/,""));}}if(parseInt(H.crossOutPriceValue)>0){N.find("#hotel_cross_sell_crossout_message").removeClass("hidden");}if(parseInt(H.percentageSavingsString)>5){N.find("#hotel_cross_sell_percent_savings_message").removeClass("hidden");}h(".photo",N).css("background-image",'url("'+s(H.urlToThumbnail)+'")');C=N.html();C=C.replace(/\{\$index\}/g,G).replace(/\{\$hotelName\}/g,decodeURIComponent(H.localizedHotelName?H.localizedHotelName:H.hotelName)).replace(/\{\$averagePrice\}/g,decodeURIComponent(H.averagePrice)).replace(/\{\$crossOutPrice\}/g,decodeURIComponent(H.crossOutPrice)).replace(/\{\$percentageSavingsString\}/g,parseInt(H.percentageSavingsString)).replace(/\{\$hotelCity\}/g,H.hotelCity).replace(/\{\$hotelProvince\}/g,(H.hotelProvince!="")?", "+H.hotelProvince:"").replace(/urlToHotelInfosite/g,decodeURIComponent(r(H.urlToHotelInfosite)));if(K){C=C.replace(/\{\$reviewRating\}/g,H.reviewRating).replace(/\{\$reviewTotal\}/g,H.reviewTotal);}h("#"+J).append(C);}}};var k=function(D){if(typeof expClientLoggingAdapter!=="undefined"){var C=D?"success":"fail";expClientLoggingAdapter.logTrxEvent(B+C);}};return{openUrl:function(C){window.open(C);},getRecommendations:function(F,C,G,H,D){if(D=="true"){var E="action=getProductRecommendations"+"&scenario=itin"+"&rule=persona-item-item-xsell"+"&useActivityData=1"+"&filterItemtla="+F+"&ruleFallback=1"+"&maxCount=8"+"&itemId="+F;}else{var E="action=getProductRecommendations"+"&scenario=checkout"+"&rule=persona-region-item-xsell"+"&maxCount=10"+"&useHotmip=true"+"&getDetails=true"+"&itemId="+F+"&itemType=1"+"&startDate="+C+"&endDate="+G+"&travellersInfo="+H+"&tripid="+w("tripid")+"&origin=CKO";}h.ajaxSetup({timeout:5000});h.ajax({url:"/Hotel-Recommender-Ajax",timeout:5000,dataType:"json",data:E,cache:false,success:function(I){if(!I.errorCode&&I.hotelModels&&I.hotelModels.length>0){h("#hotel_cross_sell_interstitial").hide();l(I.hotelModels,"hotel_cross_sell_template","hotel_cross_sell_content");h("#hotel_cross_sell_footer").show();k(true);if(D=="true"&&I.recommenderAction){s_exp_trackClick(this,"a",I.recommenderAction);}}else{h("#hotel_cross_sell_interstitial").hide();h("#hotel_cross_sell_wrapper").hide();k(false);}},error:function(I){h("#hotel_cross_sell_interstitial").hide();h("#hotel_cross_sell_wrapper").hide();k(false);}});},_getDisplayHotelRecommendation:l,submitReservation:f};})(jQuery);
/*!  generated on 2015-09-17 20:56:54.416 PDT(-0700) in 2 ms  */

/*!  served in 1 ms  */