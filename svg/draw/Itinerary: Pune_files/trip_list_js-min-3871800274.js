function deleteTrip(c,a,b){$.ajax({type:"POST",url:b,dataType:"json",data:{tripid:c,deletedTripName:encodeURI(a),parentURL:window.location.pathname},success:function(e){var d=e;if(handleRedirectRequest(d)){return;}clearMessages();if(d==null){addMessage(MESSAGE_TYPE_ERROR,MESSAGE_DELETE_FAILURE);return;}addMessages(MESSAGE_TYPE_INFO,d.messages);addMessages(MESSAGE_TYPE_ERROR,d.errors);},error:function(f,d,e){clearMessages();addMessage(MESSAGE_TYPE_ERROR,MESSAGE_DELETE_FAILURE);}});$("#menu_setting_link").focus();}function validateTripName(){var c=$.trim($("input[name='tripname']").val());var b=$("#tripname");var a=$("#tripname-error");if(c.length==0){$("#edit-name-field").addClass("edit-name-error");a.removeClass("visuallyhidden");a.removeAttr("aria-hidden");a.text(b.attr("data-validation-message"));a.offset({top:a.offset().top,left:b.offset().left});return false;}else{$("#edit-name-field").removeClass("edit-name-error");b.removeClass("edit-name-error");a.addClass("visuallyhidden");a.attr("aria-hidden","true");a.text("");return true;}}function modifyTripName(c){var b="modifyItinNameError";var a=$.trim($("input[name='tripname']").val());if(!validateTripName()){shiftFocusToTripName();return;}uitk.popups.hideInlineNotification();$("#qtip-overlay").hide();$.ajax({type:"POST",url:c,dataType:"json",data:{tripname:a,tripnote:$("#tripnote").val()},success:function(e){var d=e;if(handleRedirectRequest(d)){return;}clearMessages();if(d==null){addMessage(MESSAGE_TYPE_ERROR,b);return;}addMessages(MESSAGE_TYPE_INFO,d.messages);if(d.errors==null||d.errors.length==0){$("#title_text").text(a);$("#description_text").text($("#tripnote").val());}},error:function(f,d,e){clearMessages();addMessage(MESSAGE_TYPE_ERROR,b);}});}function submitFrequentFlyerSelection(b,g){var a=$("#frequentFlyerSelectionForm_"+b);var f=$("#frequentFlyerSelectionPopup_"+b);function d(){$(f).find(".frequentFlyerInterstitial").show();$(f).find("a.btn-close").hide();$(a).hide();}function c(){$(f).find(".frequentFlyerInterstitial").hide();$(f).find("a.btn-close").show();$(a).show();}d();$(f).find("a.btn-close").click(function(){$(a).find(".frequentFlyerError").hide();});var j=a.find("input[name=arl]").val();var h=a.find("input[name=travelerIndex]").val();var k=[];var l=a.find("input[name=numPlans]").val();for(var e=0;e<l;e++){k[e]={planID:a.find("select[name=planSelection_"+e+"] option:selected").val(),membershipNumber:a.find("input[name=membershipNumber_"+e+"]").val(),flightAirlineCode:a.find("input[name=flightAirlineCode_"+e+"]").val()};}$.ajax({url:g,data:JSON.stringify({arl:j,travelerIndex:h,plans:k}),type:"POST",contentType:"application/json",dataType:"json"}).success(function(i){if(i.responseData.status=="SUCCESS"){var n="ffMessage";var m=window.location.search;if(m.indexOf(n)==-1){if(window.location.href.indexOf("?")<0){window.location.href=window.location.href+"?ffMessage=true";}else{window.location.href=window.location.href+"&ffMessage=true";}}else{location.reload();}window.location.hash="";}else{a.find(".frequentFlyerError").show();c();}}).error(function(i){switch(i.status){case 401:window.location="${globalUrls.getLoginUrl()}&selc=0&uurl="+encodeURI(window.location);break;default:a.find(".frequentFlyerError").show();c();break;}});}$(document).ready(function(){var a=$("textarea[name='tripnote']");a.bind("keyup keydown focus blur",function(){limitLength($(this),255);});a.live("input paste",function(){limitLength($(this),255);});});function shiftFocusToTripName(){setTimeout(function(){$("#tripname").focus();},250);}function openEditTripName(){var a=$("input[name='tripname']");a.val($.trim($("#title_text").text()));validateTripName();$('a[href="#edit_modal_popup"]').click();shiftFocusToTripName();$(".btn-close").click(function(b){setTimeout(function(){$("#menu_setting_link").focus();},250);});$("#closeEditPopup").click(function(b){setTimeout(function(){$("#menu_setting_link").focus();},250);});$("#cancelEditPopup").click(function(b){setTimeout(function(){$("#menu_setting_link").focus();},250);});}$(document).ready(function(){$("span.visually-hidden").css("display","inline-block");$(".btn-close").attr("aria-label","Close");$(".frequentFlyerSelection").parent().find(".btn-close").click(a);function a(){var b=$(this).parent().find(".frequentFlyerSelection").data("control");setTimeout(function(){$("#"+b).focus();},250);}});
/*!  generated on 2015-09-18 00:11:04.279 PDT(-0700) in 3 ms  */

/*!  served in 1 ms  */