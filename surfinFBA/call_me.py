##########	Function that can call you when your script is done running.

##
##	Requires a gmail account with an app-password
##


import smtplib
carriers = {
	'att':    '@txt.att.net',
	'tmobile':' @tmomail.net',
	'verizon':  '@vtext.com',
	'sprint':   '@page.nextel.com'
}

def send(message,number,auth):
        # number is a phone number no in format XXX-XXX-XXXX
		#auth is tuple ("email@address.com", "password")

	number2 = str(number).replace('-','') + '{}'

	to_number = number2.format(carriers['att'])

	# Establish a secure session with gmail's outgoing SMTP server using your gmail account
	server = smtplib.SMTP( "smtp.gmail.com", 587 )
	server.starttls()
	server.login(auth[0], auth[1])

	# Send text message through SMS gateway of destination number
	server.sendmail( auth[0], to_number, message)
