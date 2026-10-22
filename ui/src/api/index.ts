import axios from "axios";

export type ResultData =
  | Record<string, Record<string, string | number>>
  | { data: string };

const baseURL = "/api";

const axiosInstance = axios.create({ baseURL });

const api = {
  data: {
    primerDriver(body: object) {
      return axiosInstance.post<ResultData>("", body);
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
