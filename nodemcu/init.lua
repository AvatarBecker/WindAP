-- load credentials, 'SSID' and 'PASSWORD' declared and initialize in there
dofile("credentials.lua")

function startup()
    if file.open("init.lua") == nil then
        print("init.lua deleted or renamed")
    else
        print("Running")
        file.close("init.lua")
        -- the actual application is stored in 'application.lua'
        dofile("application.lua")
    end
end

wifi.setmode(wifi.SOFTAP)
wifi.ap.config({ssid=SSID, pwd=PASSWORD, auth=wifi.WPA2_PSK})
startup()
