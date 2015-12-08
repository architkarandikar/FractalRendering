(function(){var a=50;window.AlertsMenu={_alerts:[],_endpoints:{userAlertsService:"/news/v1",expwebDelete:"/api/userhistory/clearNews"},omLock:true,queuedOmPixels:[],init:function(){var c=this,b=this.menuView.init();setTimeout(function(){c.omLock=false;c.sendQueuedOmPixels();},1000);this._getAlerts().then(function(d){AlertsMenu._unpackRawAlertData(d);AlertsMenu._translateBetaMessages();b.reRender();});this.wotifCheck().detectIE11();return this;},getAlerts:function(){return this._alerts;},count:function(){return this._alerts.length;},at:function(b){b=b||0;return this._alerts[b];},wotifCheck:function(){if(Scratchpad.state.isWotif){$(".alerts-display-wrapper").addClass("isWotif");}return this;},detectIE11:function(){if(/Trident.*rv[ :]*11\./.test(navigator.userAgent)){$(".alert-menu-wrapper").addClass("ie11");}return this;},addLateAlert:function(b){if(this.count()===0){this._alerts.push(new UserAlert({newsFeedType:b.newsFeedType||"Bad type",_uid:b._uid||++a,msgid:b.msgid||"error",params:b.params||[]}));}this.menuView.reRender();return this._alerts;},pluckAlertByUniqueId:function(d){var c=0,b=this._alerts.length;for(c;c<b;c++){if(this._alerts[c].uniqueId===d){return this._alerts.splice(c,1)[0];}}return false;},sendeVar28:function(b){if(this.omLock){this.queuedOmPixels.push(b);}else{if(window.dctk&&dctk.onReady){dctk.onReady(function(){dctk.trackImpression("APE:"+b);});}else{if(window.s_exp&&window.s_exp.tl){s_exp.eVar28="APE:"+b;s_exp.linkTrackVars="eVar28,events";s_exp.linkTrackEvents="event80,event81";s_exp.tl(true,"o","RFRR Action Link");}}}return this;},sendQueuedOmPixels:function(){for(var b=0;b<this.queuedOmPixels.length;b++){this.sendeVar28(this.queuedOmPixels[b]);}this.queuedOmPixels=[];return this;},setLocalFlag:function(b,c){try{if(window.localStorage){try{window.localStorage.setItem(b,c);}catch(d){return false;}}}catch(d){return false;}return true;},getLocalFlag:function(b){var d=false;try{if(window.localStorage){try{d=window.localStorage.getItem(b);}catch(c){return false;}}}catch(c){return false;}return d;},deleteSingle:function(c,b){if(!c){}else{this.sendeVar28("Header.Notif.Delete."+c.raw.msgid);if(b){this._remoteDelete(c);}c=null;}return this;},_getAlerts:function(){return this._request(this._buildAlertQuery());},_buildAlertQuery:function(){var e=$("#newsFeed-scalatra-app-url").attr("data-id"),g=window.alertMenuRequestIds,f=g.TUID||-1,b=$("#scratchpad-badge-guid").attr("data-id")||-1,d=g.EXPID||-1,c=g.SITEID;return e+this._endpoints.userAlertsService+"/expid/"+d+"/siteid/"+c+"/tuid/"+f+"/guid/"+b+"?format=jsonp"+"&callback=?";},_request:function(b){return $.ajax({url:b,dataType:"jsonp",complete:function(){}});},_unpackRawAlertData:function(d){if(!d.responsestatus||!d.responsestatus.status||d.responsestatus.status!=="SUCCESS"||!d.newsItems||!d.newsItems.length){return[];}this._alerts=[];var c=0,b=d.newsItems.length,e;for(c;c<b;c++){e=new UserAlert(d.newsItems[c]);if(typeof e.contentString!="undefined"){this._alerts.push(e);}}return this._alerts;},_translateBetaMessages:function(){var c=0,b=this._alerts.length;for(c;c<b;c++){this._alerts[c]._translateSPBetaMessage();}return this;},_remoteDelete:function(b){return $.ajax({url:this._endpoints.expwebDelete,type:"POST",dataType:"json",contentType:"application/json",data:JSON.stringify({type:b.raw.type||"",msgid:b.raw.msgid||""})});}};window.AlertsMenu.menuView={state:{legacyMode:false},$els:{headerBadge:null,menuWrapper:null,scrollable:null},init:function(b){b=b||$("#user-alerts-menu");this.$els.menuWrapper=b;this.state.legacyMode=!(window.uitk&&uitk.version&&uitk.version.split(".")[0]==="v2");if(b.hasClass("scratchpad-greeting")){this.prepScratchpadTray(b);}else{this.becomeMenu();}return this;},reRender:function(){this.empty();var c="",e=AlertsMenu.getAlerts(),d=0,b=AlertsMenu.count();for(d;d<b;d++){c+=e[d].toDisplayString();}this.fill(c);this.updateBadge(b);this._bindNewAlertLinks(this.state.legacyMode);this._bindNewAlertCloseButtons(this.state.legacyMode);return this;},renderToScratchpadTray:function(e){if(!e.length||e.length<1){return this;}var c="",d=0,b=e.length;for(d;d<b;d++){c+=e[d].toDisplayString();}if(c!==""){this.$els.scrollable.empty().html("<ul>"+c+"</ul>");$(".alert-item .btn-action:not(.bound)").addClass("bound").click(function(f){f.preventDefault();f.stopPropagation();Scratchpad.sendeVar28("Header.TrayAlert.SP");window.location.href="/scratchpad?rfrr=SP.Tray.LiveAlerts";});}this._bindTempNewAlertButtons(this.state.legacyMode);return this;},toggle:function(){if(this.$els.menuWrapper.find(".menu").parent().hasClass("open")){this.close();}else{this.open();}return this;},open:function(c){var b=(this.state.legacyMode)?"0.5em":"-171.5px";this.$els.menuWrapper.find(".menu").css({"right":b}).parent().addClass("open");if(c){AlertsMenu.sendeVar28("Header.Notif.MenuAuto");}else{AlertsMenu.sendeVar28("Header.Notif.MenuOpen");}return this;},close:function(){var b=(this.state.legacyMode)?"0.5em":"-195.5px";this.$els.menuWrapper.find(".menu").css({"right":b}).parent().removeClass("open");this.$els.menuWrapper.find(".bell-trigger-link").attr("aria-expanded","false");this.$els.menuWrapper.focus();AlertsMenu.sendeVar28("Header.Notif.MenuClose");return this;},empty:function(){this.$els.scrollable.find(".alert-item:not(.default-alert-item)").remove();this.$els.scrollable.find(".default-alert-item").show();return this;},fill:function(d){var b=this.$els.scrollable.find(".default-alert-item"),e=0,c=AlertsMenu.count(),f={};$(d).insertBefore(b);for(e;e<c;e++){f=AlertsMenu.at(e);f.$el=$("#"+f.uniqueId);}return this;},prepScratchpadTray:function(d){var c=Scratchpad.headerTray,b=$("<ul>").addClass("scrollable"),e=(d).empty().append(b);this.$els.headerBadge=c.$elements.$badge;this.$els.scrollable=e.find(".scrollable");c.state.greeting=true;if(c.state.currentCount!==0){$(".user-history-tab .menu").addClass("two-column");}c.$elements.$greeting.show();$(".scratchpad-item.first").removeClass("first");return this;},becomeMenu:function(){var b=this;if(this.state.legacyMode){this.$els.menuWrapper.addClass("legacy");}$(".alerts-display-wrapper").show();this.$els.scrollable=this.$els.menuWrapper.find(".scrollable");this.$els.headerBadge=this.$els.menuWrapper.find(".badge").click(function(c){c.preventDefault();c.stopPropagation();$(".shop-nav > li").removeClass("open");b.toggle();});this.$els.menuWrapper.addClass("visible").click(function(){if($(this).find(".menu").parent().hasClass("open")){AlertsMenu.sendeVar28("Header.Notif.MenuClose");}else{AlertsMenu.sendeVar28("Header.Notif.MenuOpen");}});this.updateBadge(AlertsMenu.count());$(".alerts-display-wrapper").show();this._bindNewAlertLinks(this.state.legacyMode);return this;},updateBadge:function(b){if(!b&&b!==0){b=AlertsMenu.count();}if(b>0){this.$els.headerBadge.html(b).show();if(b==1){$("#accessibleBadge").html(accessibleBadgeMessage).attr("aria-hidden","false");}else{$("#accessibleBadge").html(accessibleBadgeMessagePlural.replace("####",b)).attr("aria-hidden","false");}}else{this.$els.headerBadge.hide();$("#accessibleBadge").attr("aria-hidden","true");}return this;},_bindNewAlertLinks:function(b){var c=this;if(Scratchpad.state.loggedInUser){$(".spAccessibleButtonLink").attr("role","button");$(".spAccessibleButtonLink").on("keypress",function(d){if(d.which===32){d.preventDefault();$(this).trigger("click");}});}this.$els.menuWrapper.find(".alert-message a:not(.bound)").addClass("bound").click(function(f){var d=$(this),e=window.alertMenuClickHandlers[d.attr("data-action")],h=d.closest(".alert-item").attr("id"),g=AlertsMenu.pluckAlertByUniqueId(h);AlertsMenu.sendeVar28("Header.Notif.CTAClick");if($.isFunction(e)){e(f,g);}if(d.attr("data-trigger-delete")!=="false"){g.deleteSelf(true);c.close().updateBadge();}});return this;},_bindNewAlertCloseButtons:function(b){var c=this,d=(b)?"x":"";this.$els.menuWrapper.find(".close-button button:not(.bound)").addClass("bound").click(function(f){f.preventDefault();f.stopPropagation();var h=$(this).closest(".alert-item").attr("id"),g=AlertsMenu.pluckAlertByUniqueId(h);
AlertsMenu.sendeVar28("Header.Notif.Dismiss."+g.raw.msgid);g.deleteSelf(true);c.close().updateBadge();});return this;},_bindTempNewAlertButtons:function(b){var c=this;this.$els.menuWrapper.find(".btn-close:not(.bound)").addClass("bound").click(function(d){d.preventDefault();d.stopPropagation();var g=$(this).closest(".alert-item").attr("id"),f=AlertsMenu.pluckAlertByUniqueId(g);AlertsMenu.sendeVar28("Header.Notif.Dismiss."+f.raw.msgid);f.deleteSelf(true);c.close().updateBadge();});this.$els.menuWrapper.find(".alert-message a:not(.bound)").addClass("bound").click(function(f){f.preventDefault();f.stopPropagation();var d=$(this),e=window.alertMenuClickHandlers[d.attr("data-action")],h=d.closest(".alert-item").attr("id"),g=AlertsMenu.pluckAlertByUniqueId(h);AlertsMenu.sendeVar28("Header.Notif.SP");if($.isFunction(e)){e(f,g);}if(d.attr("data-trigger-delete")!=="false"){g.deleteSelf(true);c.close().updateBadge();}});return this;}};window.UserAlert=function(b){this.raw=b;this.contentString=this._buildContentString();this.uniqueId=[this.raw.msgid,++a].join("");this.$el=$(this.toDisplayString);return this;};UserAlert.prototype={deleteSelf:function(b){this.$el.remove();if(this.raw.newsFeedType==="scratchpadLegacyAlert"){if(this.raw.msgid=="spNudgeEnroll"){Scratchpad.declineNudge();}else{if(this.raw.msgid=="spSignIn"){Scratchpad.declineSignInToSave();}}AlertsMenu.deleteSingle(this,false);}else{AlertsMenu.deleteSingle(this,b);}return null;},toDisplayString:function(){var b="<li id="+this.uniqueId+' class="alert-item">'+this.contentString+'<div class="close-button">'+'<button id="'+this.uniqueId+'-delete" btnType="close" targetId="'+this.uniqueId+'"><span class="icon icon-close" aria-hidden="true">'+'</span><span class="visuallyhidden">'+accessibleNotifClose+"</span></button>"+"</div>"+"</li>";return b;},_buildContentString:function(){return this._fillInDynamicData(this._findContentString());},_findContentString:function(){contentString=alertMenuDisplayStrings[this.raw.msgid];return contentString;},_fillInDynamicData:function(c){if(typeof c!="undefined"){for(var b=0;b<this.raw.params.length;b++){c=c.replace("####",this.raw.params[b]);}}return c;},_translateSPBetaMessage:function(){if($.inArray(this.raw.msgid,["SPHC","SPHCONEUP","SPHCONEDOWN","SPHCBOTHONE"])===-1){return;}var b=this.raw.params[0],h=this.raw.params[1],c=(h>0),i=(h===1),e=(b>0),f=(b===1),d=(c&&e),j=(i&&f),g={msgid:"",params:[]};if(d){g.params=[b,h];if(j){g.msgid="SPHCBOTHONEV2";}else{if(f){g.msgid="SPHCONEDOWNV2";}else{if(i){g.msgid="SPHCONEUPV2";}else{g.msgid="SPHCV2";}}}}else{if(c){g.params=[h];if(i){g.msgid="SPHCONEUPONLYV2";}else{g.msgid="SPHCUPONLYV2";}}else{if(e){g.params=[b];if(f){g.msgid="SPHCONEDOWNONLYV2";}else{g.msgid="SPHCDOWNONLYV2";}}}}this.raw=g;this.contentString=this._buildContentString();return this;}};window.AlertsMenu.tests={init:function(){this.liveNudgeAlertTest.bucket=(parseInt(window.scratchpadUseLiveAlertsAsNews,10)||0);},liveNudgeAlertTest:{bucket:0,run:function(){var b=[this.control,this.variantOne][this.bucket]();return this.bucket;},control:function(){Scratchpad.headerTray._showLegacyNewsHeader();},variantOne:function(){AlertsMenu.init();if(AlersMenu.count()===0){AlertsMenu.addLateAlert({newsFeedType:"scratchpadLegacyAlert",msgid:"spBasicNews",params:[]});}}}};})();
/*!  generated on 2015-09-17 23:54:12.540 PDT(-0700) in 2 ms  */

/*!  served in 1 ms  */