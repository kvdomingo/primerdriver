import axios from "axios";

const baseURL = "/api";

const axiosInstance = axios.create({ baseURL });

const api = {
  data: {
    primerDriver(body: object) {
      return axiosInstance.post("", body);
    },
    version() {
      return axiosInstance.get<{ program_version: string; web_version: string }>(
        "/version",
      );
    },
    expressionSystem() {
      return axiosInstance.get<{ data: string[] }>("/expressionsys");
    },
  },
};

export default api;
