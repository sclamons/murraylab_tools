import smtplib
from email.message import EmailMessage
import warnings
import __main__ as main
import os
import traceback
import socket

# We're going to use some GLOBAL VARIABLES to keep track of user information.
# This should be set once somewhere at the beginning of the program
# using the set_credentials method below.
_username = None
_password = None
_destination = None

def set_credentials(username, password):
    if not isinstance(username, str):
        raise TypeException("Username must be a string.")
    if not isinstance(password, str):
        raise TypeException("Password must be a string.")
    global _username
    global _password
    _username = username
    _password = password

def set_destination(destination):
    if not isinstance(destination, str):
        raise TypeException("Destination must be a string.")
    global _destination
    _destination = destination

class notify_when_done(object):
    def __init__(self, done_string, debug_connection = False,
                  server = "smtp.gmail.com", port = 587):
        self.done_string = done_string
        self.debug_connection = debug_connection
        self.server = server
        self.port   = port

        if hasattr(main, "__file__"):
            self.file = main.__file__
        else:
            self.file = "Jupter notebook at " + os.path.abspath("")

    def __enter__(self):
        # Test connection and credentials.
        self.login_good = True
        try:
            server = self.connect_to_server()
        except Exception as e:
            warnings.warn("Connection to email server failed. "\
                          "Run will continue without notifications.",
                          RuntimeWarning)
            self.login_good = False
            if self.debug_connection:
                print("Connection error is:")
                print(traceback.format_exc())

    def __exit__(self, type, value, trace):
        # Don't bother connecting again if it didn't work the first time.
        if not self.login_good:
            return

        # Send an email/text informing the owner of the script what happened.
        server = self.connect_to_server()

        message = self.file + " reports "
        if type == None:
            # Block completed without error.
            message += "completion of: " + self.done_string
        else:
            # Block hit an error.
            message += "the following error: \n" + traceback.format_exc()

        msg_email = EmailMessage()
        msg_email.set_content(message)
        msg_email["Subject"] = ("Error " if type != None else "Success ") +\
                                "reported by " + self.file
        msg_email["From"]    = _username
        msg_email["To"]      = _destination

        server.send_message(msg_email)

    def connect_to_server(self):
        # Use socket.gethostbyname to get an IPv4 address, which connects much,
        # much faster than an IPv6 address.
        server = smtplib.SMTP(socket.gethostbyname(self.server), self.port)
        server.starttls()

        server.login(_username, _password)
        return server