import gspread
from oauth2client.service_account import ServiceAccountCredentials

# Authenticate with Google Sheets API using credentials
scope = ["https://spreadsheets.google.com/feeds", "https://www.googleapis.com/auth/drive"]
creds = ServiceAccountCredentials.from_json_keyfile_name("your-credentials.json", scope)
client = gspread.authorize(creds)

# Open the Google Sheet by its title or URL
sheet = client.open("Your Google Sheet Name").sheet1

# Get the data from the sheet
data = sheet.get_all_values()

# Save the data as a CSV file
with open("output.csv", "w") as csv_file:
    for row in data:
        csv_file.write(",".join(row) + "\n")