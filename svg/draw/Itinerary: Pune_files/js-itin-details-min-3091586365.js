var TripDetails=(function(){function b(c){c=c||{};this.dataTrackAttribute=c.dataTrackAttribute||"data-track";this.interstitialLinkClass="."+(c.interstitialLinkClass||"td_interstitial_link");this.populateSpinnerForm=c.populateSpinnerForm||populateSpinnerForm;}b.prototype={init:function(){$(this.interstitialLinkClass).click(a(this));}};function a(c){return function(){var d=$(this).attr("href"),e=$(this).attr(c.dataTrackAttribute);c.populateSpinnerForm(d);uitk.popups.displayInterstitialModal("#spinner","#spinnerForm");s_exp_trackClick(this,"a",e);return false;};}return b;})();
/*!  generated on 2015-09-18 00:16:09.807 PDT(-0700) in 1 ms  */

/*!  served in 0 ms  */
