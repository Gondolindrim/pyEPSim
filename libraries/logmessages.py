import colorama
from colorama import Fore, Back, Style
colorama.init()

# ERROR CHECK FUNCTIONS AND ERRORS ------------------------------------------------------------ {{{1
# These functions and exceptions are custom-tailored to handle error checks in the loadCase function.
# error_message is a string function that formats the string in a constrasted format to indicate an error has ocurred: bold, red "ERROR" with an arrow and a string.
# warning_message and good_message have the same spirit, but in yellow and green.
def error_message(string): return Style.BRIGHT + Fore.RED + ' -> ERROR: ' + Style.RESET_ALL + string
def warning_message(string): return '\n'+ Style.BRIGHT + Fore.YELLOW + ' -> WARNING: ' + Style.RESET_ALL + string
def good_message(string): return Style.BRIGHT + Fore.GREEN + ' -> ' + Style.RESET_ALL + string
def header_message(string): return Style.BRIGHT + Fore.WHITE + 50*'-' + '\n' + string + '\n' + 50*'-' + Style.RESET_ALL
def green_string(string): return Style.BRIGHT + Fore.GREEN + string + Style.RESET_ALL
def green_indicator_string(indicator, string): return Style.BRIGHT + Fore.GREEN + indicator + Fore.WHITE + string + Style.RESET_ALL
def bold_string(string): return Style.BRIGHT + string + Style.RESET_ALL

