   var inauth = inauth || {};
    inauth.getCookie = function (c_name) {
        var i, x, y, ARRcookies = document.cookie.split(";");
        for (i = 0; i < ARRcookies.length; i++) {
            x = ARRcookies[i].substr(0, ARRcookies[i].indexOf("="));
            y = ARRcookies[i].substr(ARRcookies[i].indexOf("=") + 1);
            x = x.replace(/^\s+|\s+$/g, "");
            if (x == c_name) {
                return unescape(y);
            }
        }
    }
    inauth.setCookie = function (c_name, value, exdays) {
        var exdate = new Date();
        exdate.setDate(exdate.getDate() + exdays);
        var c_value = escape(value) + ((exdays == null) ? "" : "; expires=" + exdate.toUTCString());
        document.cookie = c_name + "=" + c_value;
    }

    inauth.init = function () {
        var mc1 = inauth.getCookie('MC1');
        if (mc1 && mc1.length > 5 && mc1.indexOf('GUID=') == 0) {
            inauth.guid = mc1.substr(5);
        }
        inauth.guid = inauth.guid || "";

        inauth.eid = inauth.getCookie('eid') || "";
        inauth.IAID = inauth.getCookie('IAID') || "";
    }
    inauth.init();

if(!inauth.IAID){
    var _cc = _cc || [];
_cc.push(['st', 1000]);
    _cc.push(['ci', { 'sid': 'eeb0e3596d8e435b89bf23a76', 'tid': Math.random().toString(36).substr(2, 9), 'expuserid': inauth.eid, 'guid': inauth.guid }]);
    _cc.push(['crdi', function (deviceInfo) {
        var print = deviceInfo.dp;
        var resilientGuid = deviceInfo.drg;
        // store the resilientGuid in a cookie named IAID
        inauth.setCookie('IAID', resilientGuid);
    }]);

    _cc.push(['run', ('https:' == document.location.protocol ? 'https://' : 'http://') + 'www.cdn-net.com']);

    (function () {
        var c = document.createElement('script'); c.type = 'text/javascript'; c.async = true;
        c.src = ('https:' == document.location.protocol ? 'https://' : 'http://') + 'www.cdn-net.com/cc.js';
        var s = document.getElementsByTagName('script')[0]; s.parentNode.insertBefore(c, s);
    })();
}