using RydbergEmulator: Web

Web.WEBCONFIG[:PORT] = 8080
SERVER = run_server(Web.naive_handler)
