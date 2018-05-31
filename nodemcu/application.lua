-- Receives duty cycle (from 0 to 1023) per UDP datagram (one value per packet. Yeah, I know, judge me)
-- Note: to send a UDP package from command line in linux, do:
-- sudo sendip -v -p ipv4 -is 192.168.4.2 -p udp -us 3031 -ud 3031 192.168.4.1 -d 512
-- the last argument being the duty cycle

local led = 4
local duty_cycle = 1023
local udp_socket = net.createUDPSocket()

-- Routine to update PWM duty cycle
local function DutyUpdate(source_socket, duty, port, ip)
    duty_cycle = duty
    pwm.setduty(led, duty_cycle)
    print("New duty cycle: "..duty_cycle)
    print("Data received from "..ip)
end

-- Define WiFi station event callbacks 
WifiConnectEvent = function(T) 
    pwm.setduty(led, duty_cycle)
    print("Connection to STA("..T.MAC..") established!")
    print("LED ON! Resuming previous duty cycle")
end

WifiDisconnectEvent = function(T)
    pwm.setduty(led, 1023)
    print("STA("..T.MAC..") disconnected!")
    print("LED OFF!")
end
-- pwm at pin "led", 1 kHz (shame on you NodeMCU!), with 100% duty cycle
pwm.setup(led, 1000, duty_cycle)
pwm.start(4)

-- Register WiFi Station event callbacks (wifi stuff created in init)
wifi.eventmon.register(wifi.eventmon.AP_STACONNECTED, WifiConnectEvent)
wifi.eventmon.register(wifi.eventmon.AP_STADISCONNECTED, WifiDisconnectEvent)

udp_socket:listen(3031, wifi.ap.getip())

-- each time a UDP datagram is received, update the PWM duty cycle
udp_socket:on("receive", DutyUpdate)
